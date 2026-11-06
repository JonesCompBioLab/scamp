from sklearn.mixture import GaussianMixture
from kneed import KneeLocator
import numpy as np
import torch
import math
from scipy.stats import norm
import hashlib


def _vector_seed(input):
    """Random seed derived from a copy-number vector, so reruns are identical.

    The seed is the first 8 bytes of a sha256 of the values rounded to 3 decimals
    (the rounding `predict` uses to spot duplicate vectors; + 0.0 turns -0.0 into 0.0).
    """
    rounded = np.round(np.asarray(input, dtype=np.float64).ravel(), 3) + 0.0
    digest = hashlib.sha256(rounded.tobytes()).digest()
    return int.from_bytes(digest[:8], "little")


def _trim_cells(input, max_percentile=100.0, min_copy_number=None):
    """Keeps the cells at or below `max_percentile` of the copy-numbers and strictly
    above `min_copy_number` (None = no minimum). Both are applied to all of the
    gene's cells. Returns a 1-D array.
    """
    x = np.asarray(input, dtype=np.float64).ravel()
    keep = x <= np.percentile(x, max_percentile)
    if min_copy_number is not None:
        keep &= x > min_copy_number
    return x[keep]


def GMM_ecDNA_predict(
    model,
    input,
    out_log,
    n_components=None,
    one_thresh=1.15,
    max_components=15,
    genes=np.array(["GENE"]),
    decision_rule=0.5,
    kneedle_coeff=3,
    var_scale=0.5,
    min_weight=0.005,
    verbose=False,
    gmm_max_percentile=100.0,
    gmm_min_copy_number=None,
    nn_min_copy_number=2.0,
    nn_max_percentile=99.0,
    nn_filter_copy_number=2.5,
    score="ratio",
    n_draws=1,
    single_dist_scamp=False,
    scamp_min_copy_number=2.0,
    scamp_max_percentile=99.9,
    scamp_filter_copy_number=2.5,
):
    """Predicts ecDNA status for one gene by fitting a GMM to its copy-numbers.

    Each fitted distribution with enough weight is sampled (with a seed derived
    from the copy-numbers and the distribution's index), and the samples are
    classified with `NN_ecDNA_predict`. The gene is ecDNA-positive if any
    distribution is, and the first positive distribution is returned.

    Args:
        model: SCAMP model loaded from a .pt file.
        input: Copy-numbers of one gene across cells.
        out_log: List of log lines; messages are appended to it when `verbose`.
        n_components: Number of GMM distributions, if known. If None, it is
            chosen with kneedle.
        one_thresh: If `n_components` is None, a single distribution is used when
            the 1-distribution BIC divided by the minimum BIC is below this.
        max_components: If `n_components` is None, fits from 1 up to
            `max_components - 1` distributions.
        genes: Array holding the gene name, used in log messages.
        decision_rule: SCAMP likelihood threshold for calling a distribution
            ecDNA-positive.
        kneedle_coeff: Higher values make kneedle pick fewer distributions.
        var_scale: Distribution variances are scaled as var * weight ** var_scale
            before sampling, to counteract small sample sizes.
        min_weight: Distributions with less weight than this (out of 1) are ignored.
        verbose: Whether to log per-distribution details.
        gmm_max_percentile: Before the GMM fit, drop cells above this percentile
            of the gene's copy-numbers (100 = keep all).
        gmm_min_copy_number: Before the GMM fit, keep only cells with copy-number
            strictly above this (None = no minimum).
        nn_min_copy_number, nn_max_percentile, nn_filter_copy_number: Passed to
            `NN_ecDNA_predict` (the model's preprocessing of what it scores).

    Returns:
        A tuple `(prediction, probability, mean, var, weight)`: the predicted ecDNA
        status, and the SCAMP probability, mean, unscaled variance and weight of the
        first ecDNA-positive distribution. If no distribution is positive, returns
        `(False, probability, -1, -1, -1)`, where `probability` is the highest
        among the distributions checked (0 if none was).
    """

    input = np.array(input).reshape(-1, 1)
    # Seeded from the gene's copy-numbers, so the sampled cells are the same on every run
    seed = _vector_seed(input)
    all_cells = input.copy()   # untrimmed, for single_dist_scamp

    # Optional trimming of the real cells before the GMM fit
    if gmm_max_percentile < 100 or gmm_min_copy_number is not None:
        input = _trim_cells(
            input, gmm_max_percentile, gmm_min_copy_number
        ).reshape(-1, 1)
        # Too few cells left to fit up to max_components - 1 distributions
        if len(input) < max_components:
            out_log.append(
                f"{genes[0]}: only {len(input)} cells left after GMM trimming, not called"
            )
            return False, 0.0, -1, -1, -1

    gmm = _GMM_fit(
        input, n_components, one_thresh, max_components, kneedle_coeff
    )

    if verbose:
        out_log.append(f"Gene: {genes[0]}")

    # Only one distribution passes min_weight: score the gene as published scAmp does, on all of its
    # (untrimmed) cells, with scAmp's preprocessing and softmax >= decision_rule
    if single_dist_scamp and int((gmm.weights_ > min_weight).sum()) == 1:
        i = int(np.argmax(gmm.weights_))
        probability = _nn_probability(
            model, all_cells, genes, "softmax", scamp_min_copy_number,
            scamp_max_percentile, scamp_filter_copy_number)
        if verbose:
            out_log.append(f"Gene: {genes[0]}: one distribution above min_weight, scAmp scoring: {probability}")
        if probability >= decision_rule:
            return True, probability, gmm.means_[i, 0], gmm.covariances_[i, 0, 0], gmm.weights_[i]
        return False, probability, -1, -1, -1

    # Highest probability among the distributions, reported if none is positive
    max_probability = 0.0
    for i in range(gmm.n_components):
        mean = gmm.means_[i, 0]
        var = gmm.covariances_[i, 0, 0]
        weight = gmm.weights_[i]
        var_scaled = var * (weight**var_scale)

        if verbose:
            out_log.append(
                f"Mean: {mean}, Var: {var}, Var (scaled) : {var_scaled}, Weight: {weight}"
            )

        if weight > min_weight:
            # Draw 1 uses seed [seed, i]; draws 2.. use [seed, i, k]. The call uses the mean probability.
            probabilities = []
            for k in range(n_draws):
                samples = np.random.default_rng([seed, i] if k == 0 else [seed, i, k]).normal(
                    mean, np.sqrt(var_scaled), size=1000
                ).tolist()
                probabilities.append(_nn_probability(
                    model, samples, genes, score, nn_min_copy_number,
                    nn_max_percentile, nn_filter_copy_number, out_log, verbose))
            probability = float(np.mean(probabilities))
            prediction = probability >= decision_rule if score == "softmax" else probability > decision_rule
            if prediction:
                return True, probability, mean, var, weight
            max_probability = max(max_probability, probability)

    return False, max_probability, -1, -1, -1


def _GMM_fit(
    input,
    n_components=None,
    one_thresh=1.15,
    max_components=15,
    kneedle_coeff=3,
):
    if n_components is None:
        bic_scores = []
        aic_scores = []
        ks = []
        for k in range(1, max_components):
            gmm = GaussianMixture(n_components=k, random_state=0)
            gmm.fit(input)

            bic_scores.append(gmm.bic(input))
            aic_scores.append(gmm.aic(input))
            ks.append(k)

        # if bic_scores[0] / np.min(bic_scores) < one_thresh:
        #     n_components = 1
        
        # Relative gap between the 1-distribution BIC and the best BIC. Same as
        # bic_scores[0] / min - 1 when BIC is positive, but keeps its sign when BIC is negative.
        best_bic = np.min(bic_scores)
        if (bic_scores[0] - best_bic) / abs(best_bic) < one_thresh - 1:
             n_components = 1

        else:

            bic_scaled = (bic_scores - (np.min(bic_scores))) / np.ptp(
                bic_scores
            )
            bic_scaled = bic_scaled**kneedle_coeff
            # bic_scaled = np.log(np.array(bic_scores) - np.min(bic_scores) + 1e-8)

            kneedle = KneeLocator(
                ks, bic_scaled, curve="convex", direction="decreasing"
            )

            n_components = kneedle.knee

            # kneedle returns None when the BIC curve has no knee
            if n_components is None:
                n_components = 1

    # Fit GMM
    gmm = GaussianMixture(n_components=n_components, random_state=0)
    gmm.fit(input)

    return gmm


def _nn_probability(model, cells, genes, score, min_copy_number, max_percentile,
                    filter_copy_number, out_log=None, verbose=False):
    """The model's ecDNA probability for one set of cells: "ratio" = scaup's NN_ecDNA_predict score,
    "softmax" = published scAmp's model.proba. 0.0 if the cells fail the model's copy-number filter."""
    if score == "ratio":
        return float(NN_ecDNA_predict(
            model, cells, [] if out_log is None else out_log, genes, 0.5, verbose,
            min_copy_number=min_copy_number, max_percentile=max_percentile,
            filter_copy_number=filter_copy_number)[1])
    X, genes_pass_filter = model.prepare_copy_numbers(
        np.array(cells, dtype=float).reshape(-1, 1), genes, min_copy_number=min_copy_number,
        max_percentile=max_percentile, filter_copy_number=filter_copy_number)
    if len(genes_pass_filter) == 0:
        return 0.0
    with torch.no_grad():
        return float(model.proba(torch.Tensor(X))[0, 1])


def NN_ecDNA_predict(
    model,
    input,
    out_log,
    genes=np.array(["GENE"]),
    decision_rule=0.5,
    verbose=False,
    min_copy_number=2.0,
    max_percentile=99.0,
    filter_copy_number=2.5,
):
    """Predicts ecDNA status for one copy-number distribution with the SCAMP model.

    Args:
        model: SCAMP model loaded from a .pt file.
        input: Copy-numbers of one gene across cells.
        out_log: List of log lines; messages are appended to it when `verbose`.
        genes: Array holding the gene name, used in log messages.
        decision_rule: SCAMP likelihood threshold for calling the distribution
            ecDNA-positive.
        verbose: Whether to log when the distribution fails the copy-number filter.
        min_copy_number: The model's summary statistics ignore copy-numbers at or
            below this (-inf = keep all).
        max_percentile: Copy-numbers above this percentile are dropped before the
            summary statistics (100 = keep all).
        filter_copy_number: The distribution is not scored (False, 0) if its mean
            copy-number is below this (-inf = always score).

    Returns:
        A tuple `(prediction, probability)`: the predicted ecDNA status, and the
        model's ecDNA probability as a scalar tensor. If the distribution fails the
        model's copy-number filter, returns `(False, 0)`.
    """
    input = np.array(input).reshape(-1, 1)

    X, genes_pass_filter = model.prepare_copy_numbers(
        input,
        genes,
        min_copy_number=min_copy_number,
        max_percentile=max_percentile,
        filter_copy_number=filter_copy_number,
    )

    if len(genes_pass_filter) == 0:
        if verbose:
            out_log.append(f" {genes[0]} did not pass filter")
        return False, 0

    # Just get the second to last layer
    _x = model.forward(torch.Tensor(X))
    _x = _x[0]
    ecDNA_pred = _x[1] / (_x[0] + _x[1])

    return (ecDNA_pred > decision_rule).item(), ecDNA_pred


def KNN_ecDNA_predict(
    model,
    input,
    out_log,
    genes=np.array(["GENE"]),
    n_components=None,
    one_thresh=1.15,
    max_components=15,
    kneedle_coeff=3,
    k=100,
    ecDNA_percentage_thresh=0.001,
    var_scale=0.5,
    decision_rule=0.5,
    verbose=False,
    gmm_max_percentile=100.0,
    gmm_min_copy_number=None,
    nn_min_copy_number=2.0,
    nn_max_percentile=99.0,
    nn_filter_copy_number=2.5,
    min_weight=0.005,
    knn_min_mean=2.5,
    knn_min_var=10,
    score="ratio",
    single_dist_scamp=False,
    scamp_min_copy_number=2.0,
    scamp_max_percentile=99.9,
    scamp_filter_copy_number=2.5,
):
    """Predicts ecDNA status for one gene from the neighborhoods of individual cells.

    Cells are assigned to distributions of a GMM fit, and each cell's neighbors in
    the same distribution are classified with `NN_ecDNA_predict`. Inputs with more
    than 4000 cells are randomly subsampled to 4000, with a seed derived from the
    copy-numbers.

    Args:
        model: SCAMP model loaded from a .pt file.
        input: Copy-numbers of one gene across cells.
        out_log: List of log lines; messages are appended to it when `verbose`.
        genes: Array holding the gene name, used in log messages.
        n_components: Number of GMM distributions, if known. If None, it is
            chosen with kneedle.
        one_thresh: If `n_components` is None, a single distribution is used when
            the 1-distribution BIC divided by the minimum BIC is below this.
        max_components: If `n_components` is None, fits from 1 up to
            `max_components - 1` distributions.
        kneedle_coeff: Higher values make kneedle pick fewer distributions.
        k: Number of neighbors per cell.
        ecDNA_percentage_thresh: The gene is ecDNA-positive if more than this
            fraction of cells are ecDNA-positive.
        var_scale: Distribution variances are scaled as var * weight ** var_scale
            when assigning cells, to counteract small sample sizes.
        decision_rule: SCAMP likelihood threshold for calling a neighborhood
            ecDNA-positive.
        verbose: Whether to log per-cell details.
        gmm_max_percentile: Drop cells above this percentile of the gene's
            copy-numbers (100 = keep all). The trimmed cells are used for the GMM
            fit, the subsample and the neighborhoods.
        gmm_min_copy_number: Keep only cells with copy-number strictly above this
            (None = no minimum). Used like `gmm_max_percentile`.
        nn_min_copy_number, nn_max_percentile, nn_filter_copy_number: Passed to
            `NN_ecDNA_predict` (the model's preprocessing of what it scores).

    Returns:
        A tuple `(prediction, num_positive_cells)`: the predicted ecDNA status, and
        the number of cells whose neighborhood was called ecDNA-positive.
    """

    # Seeded from the gene's copy-numbers, so the subsample is the same on every run
    seed = _vector_seed(input)
    all_cells = np.asarray(input, dtype=float).copy()   # untrimmed, for single_dist_scamp

    # Optional trimming of the real cells, before the subsample and the GMM fit
    if gmm_max_percentile < 100 or gmm_min_copy_number is not None:
        input = _trim_cells(input, gmm_max_percentile, gmm_min_copy_number)
        # Too few cells left to fit up to max_components - 1 distributions
        if len(input) < max_components:
            out_log.append(
                f"{genes[0]}: only {len(input)} cells left after GMM trimming, not called"
            )
            return False, 0

    if len(input) > 4000:
        arr = np.random.default_rng(seed).choice(input, size=4000, replace=False)
    else:
        arr = np.asarray(input)

    # Sort the cells by copy-number so neighbours in the sliding window are the nearest
    # copy-numbers. The sort was lost in 6ae0f12; Kevin's paper code does arr = sorted(arr).
    arr = np.sort(arr)

    if verbose:
        out_log.append(f"k = {k}")

    gmm = _GMM_fit(
        np.array(input).reshape(-1, 1),
        n_components,
        one_thresh,
        max_components,
        kneedle_coeff,
    )
    x = np.array(arr).reshape(-1)

    # Only one distribution passes min_weight: score the gene as published scAmp does (P8, as in GMM
    # mode): all of its untrimmed cells, scAmp's preprocessing, softmax >= decision_rule.
    # NumberOfPositiveCells is then reported as -1 (no neighbourhoods scored).
    if single_dist_scamp and int((gmm.weights_ > min_weight).sum()) == 1:
        probability = _nn_probability(
            model, all_cells, genes, "softmax", scamp_min_copy_number,
            scamp_max_percentile, scamp_filter_copy_number)
        return probability >= decision_rule, -1

    means = gmm.means_.flatten()
    variances = gmm.covariances_.flatten()
    weights = gmm.weights_

    # Adjust variance by weight
    adjusted_variances = variances * (weights**var_scale)

    scores = []
    for mean, var, weight in zip(means, adjusted_variances, weights):
        scores.append(
            np.log(weight) + norm.logpdf(x, loc=mean, scale=np.sqrt(var))
        )

    scores = np.vstack(scores)

    # pick most likely component
    labels = np.argmax(scores, axis=0)

    n = len(arr)
    out = []

    window_size = k + 1  # include itself

    # Get label means to normalize distributions
    arr = np.asarray(arr)
    label_means = {
        label: arr[labels == label].mean() for label in np.unique(labels)
    }

    # Slide window
    for i in range(n):
        target = arr[i]
        target_label = labels[i]

        l, r = i - 1, i + 1
        count = 1

        while count < window_size:
            left_ok = l >= 0 and labels[l] == target_label
            right_ok = r < n and labels[r] == target_label

            if not left_ok and not right_ok:
                break

            if left_ok and right_ok:
                if abs(arr[l] - target) <= abs(arr[r] - target):
                    l -= 1
                else:
                    r += 1
            elif left_ok:
                l -= 1
            else:
                r += 1

            count += 1

        subset = np.array(arr[l + 1 : r])
        mean_diff = np.mean(subset) - label_means[target_label]
        subset = subset - mean_diff

        # Gate: Kevin set mean > 2.5 and var > 10 to save runtime (P8: configurable; -inf = off)
        if np.mean(subset) > knn_min_mean and np.var(subset) > knn_min_var:
            if score == "ratio":
                res = NN_ecDNA_predict(
                    model,
                    subset.tolist(),
                    out_log,
                    genes,
                    decision_rule,
                    verbose,
                    min_copy_number=nn_min_copy_number,
                    max_percentile=nn_max_percentile,
                    filter_copy_number=nn_filter_copy_number,
                )[0]
            else:   # softmax >= decision_rule, as published scAmp
                res = _nn_probability(
                    model, subset.tolist(), genes, "softmax", nn_min_copy_number,
                    nn_max_percentile, nn_filter_copy_number, out_log, verbose) >= decision_rule
            if verbose:
                out_log.append(f"{genes[0]} ecDNA:")
            out.append(res)

    return sum(out) > (ecDNA_percentage_thresh * n), sum(out)
