from sklearn.mixture import GaussianMixture
from kneed import KneeLocator
import numpy as np
from .NN_model import SCAMP
import torch
import math
from scipy.stats import norm


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
):
    """Predicts ecDNA status for one gene by fitting a GMM to its copy-numbers.

    Each fitted distribution with enough weight is sampled, and the samples are
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

    Returns:
        A tuple `(prediction, probability, mean, var, weight)`: the predicted ecDNA
        status, and the SCAMP probability, mean, unscaled variance and weight of the
        first ecDNA-positive distribution. If no distribution is positive, returns
        `(False, probability, -1, -1, -1)`, where `probability` is the highest
        among the distributions checked (0 if none was).
    """

    input = np.array(input).reshape(-1, 1)

    gmm = _GMM_fit(
        input, n_components, one_thresh, max_components, kneedle_coeff
    )

    if verbose:
        out_log.append(f"Gene: {genes[0]}")

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
            samples = np.random.normal(
                mean, np.sqrt(var_scaled), size=1000
            ).tolist()
            prediction, probability = NN_ecDNA_predict(
                model, samples, out_log, genes, decision_rule, verbose
            )
            probability = float(probability)
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


def NN_ecDNA_predict(
    model,
    input,
    out_log,
    genes=np.array(["GENE"]),
    decision_rule=0.5,
    verbose=False,
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

    Returns:
        A tuple `(prediction, probability)`: the predicted ecDNA status, and the
        model's ecDNA probability as a scalar tensor. If the distribution fails the
        model's copy-number filter, returns `(False, 0)`.
    """
    input = np.array(input).reshape(-1, 1)

    X, genes_pass_filter = model.prepare_copy_numbers(
        input,
        genes,
        min_copy_number=2.0,
        max_percentile=99.0,
        filter_copy_number=2.5,
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
):
    """Predicts ecDNA status for one gene from the neighborhoods of individual cells.

    Cells are assigned to distributions of a GMM fit, and each cell's neighbors in
    the same distribution are classified with `NN_ecDNA_predict`. Inputs with more
    than 4000 cells are randomly subsampled to 4000.

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

    Returns:
        A tuple `(prediction, num_positive_cells)`: the predicted ecDNA status, and
        the number of cells whose neighborhood was called ecDNA-positive.
    """

    if len(input) > 4000:
        arr = np.random.choice(input, size=4000, replace=False)
    else:
        arr = np.asarray(input)

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

        if np.mean(subset) > 2.5 and np.var(subset) > 10:
            res = NN_ecDNA_predict(
                model, subset.tolist(), out_log, genes, decision_rule, verbose
            )[0]
            if verbose:
                out_log.append(f"{genes[0]} ecDNA:")
            out.append(res)

    return sum(out) > (ecDNA_percentage_thresh * n), sum(out)
