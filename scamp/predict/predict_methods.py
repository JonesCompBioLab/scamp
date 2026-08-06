from sklearn.mixture import GaussianMixture
from kneed import KneeLocator
import numpy as np
from .NN_model import SCAMP
import torch
import math
from scipy.stats import norm


'''
PARAMETERS
model : SCAMP model loaded from .pt file
input : copy number list
out_log : where to send outputs to if verbose
n_components : if known, the number of components. Set to None to let kneed find
one_thresh : if n_components is None, uses this to determine when the knee is at 1 
max_components : if n_components is None, checks up to these components
genes : gene name
decision_rule : scAmp NN threshold for classifying ecDNA status
kneedle_coeff : higher means kneedle tends to smaller distributions
var_scale : changes the estimated variance to var * weight ** (var_scale), to counteract small sample size
min_weight : ignore distributions with less than this weight (max 1)

RETURNS
Boolean, the predicted ecDNA status
Mean of the true ecDNA status distribution
Var of the true ecDNA status distribution
Weight of the true ecDNA status distribution
'''
def GMM_ecDNA_predict(model, input, out_log,
                      n_components = None, 
                      one_thresh = 1.15, 
                      max_components = 15,
                      genes = np.array(['GENE']), 
                      decision_rule = 0.5,
                      kneedle_coeff = 3,
                      var_scale = 0.5,
                      min_weight = 0.005,
                      verbose = False
                      ) :
    
    input = np.array(input).reshape(-1, 1)

    gmm = _GMM_fit(input, n_components, one_thresh, max_components, kneedle_coeff)

    if verbose :
        out_log.append(f"Gene: {genes[0]}")
    
    for i in range(gmm.n_components):
        mean = gmm.means_[i, 0]
        var = gmm.covariances_[i, 0, 0]
        weight = gmm.weights_[i]
        var_scaled = var * (weight ** var_scale)

        if verbose :
            out_log.append(f"Mean: {mean}, Var: {var}, Var (scaled) : {var_scaled}, Weight: {weight}")

        if weight > min_weight :
            samples = np.random.normal(mean, np.sqrt(var_scaled), size=1000).tolist()
            if NN_ecDNA_predict(model, samples, out_log, genes, decision_rule, verbose)[0] :
                return True, mean, var, weight

    return False, -1, -1, -1


def _GMM_fit(input, n_components = None, one_thresh = 1.15, max_components = 15, kneedle_coeff = 3) :
    if n_components is None :
        bic_scores = []
        aic_scores = []
        ks = []
        for k in range(1, max_components):
            gmm = GaussianMixture(n_components=k, random_state=0)
            gmm.fit(input)

            bic_scores.append(gmm.bic(input))
            aic_scores.append(gmm.aic(input))
            ks.append(k)

        if bic_scores[0]/np.min(bic_scores) < one_thresh :
            n_components = 1

        else :

            bic_scaled = (bic_scores - (np.min(bic_scores))) / np.ptp(bic_scores)
            bic_scaled = bic_scaled ** kneedle_coeff
            # bic_scaled = np.log(np.array(bic_scores) - np.min(bic_scores) + 1e-8)
            
            kneedle = KneeLocator(
                ks,
                bic_scaled,
                curve="convex",
                direction="decreasing"
            )

            n_components = kneedle.knee

    # Fit GMM
    gmm = GaussianMixture(n_components=n_components, random_state=0)
    gmm.fit(input)

    return gmm


'''
PARAMETERS
model : SCAMP model loaded from .pt file
input : copy number list
genes : gene name
decision_rule : threshold for classifying ecDNA status
verbose : Whether or not to print results

RETURNS
Boolean, the predicted ecDNA status
Proba, the probability given by scamp
'''
def NN_ecDNA_predict(model, input, out_log, genes = np.array(['GENE']), decision_rule = 0.5, verbose = False) :
    input = np.array(input).reshape(-1, 1)

    X, genes_pass_filter = model.prepare_copy_numbers(
        input,
        genes,
        min_copy_number=2.0,
        max_percentile=99.0,
        filter_copy_number=2.5,
    )


    if len(genes_pass_filter) == 0 :
        if verbose :
            out_log.append(f" {genes[0]} did not pass filter")
        return False

    # Just get the second to last layer
    _x = model.forward(torch.Tensor(X))
    _x = _x[0]
    ecDNA_pred = _x[1]/(_x[0] + _x[1])

    return (ecDNA_pred > decision_rule).item(), ecDNA_pred


'''
PARAMETERS
model : SCAMP model loaded from .pt file
input : copy number list
out_log : output log (for NN)
genes : gene name
n_components : if known, the number of components. Set to None to let kneed find
one_thresh : if n_components is None, uses this to determine when the knee is at 1 
max_components : if n_components is None, checks up to these components
kneedle_coeff : higher means kneedle tends to fewer distributions
k_mult : multiplier for k. k = k_mult * sqrt(num cells)
ecDNA_percentage_thresh : threshold for ecDNA positive cells classifying ecDNA status
var_scale : changes the estimated variance to var * weight ** (var_scale), to counteract small sample size
decision_rule : scAmp NN threshold for classifying ecDNA status
verbose: True to add to the out log

RETURNS
Boolean, the predicted ecDNA status
''' 
def KNN_ecDNA_predict(model, input, out_log, 
                      genes = np.array(['GENE']), 
                      n_components = None, 
                      one_thresh = 1.15, 
                      max_components = 15,
                      kneedle_coeff = 3,
                      k_mult = 2.5, 
                      k = None,
                      ecDNA_percentage_thresh = 0.001, 
                      var_scale = 0.5,
                      decision_rule = 0.5,
                      verbose = False) :
    arr = sorted(input)

    if k is None :
        k = int(k_mult * math.sqrt(len(input)))
        # K should never be less than 10
        k = max(k, 10)
        # K is never greater than the input
        k = min(k, len(input))
    if verbose :
        out_log.append(f"k = {k}")

    gmm = _GMM_fit(np.array(input).reshape(-1, 1), n_components, one_thresh, max_components, kneedle_coeff)
    x = np.array(arr).reshape(-1)

    means = gmm.means_.flatten()
    variances = gmm.covariances_.flatten()
    weights = gmm.weights_

    # Adjust variance by weight
    adjusted_variances = variances * (weights ** var_scale)

    scores = []
    for mean, var, weight in zip(means, adjusted_variances, weights):
        scores.append(
            np.log(weight) +
            norm.logpdf(x, loc=mean, scale=np.sqrt(var))
        )

    scores = np.vstack(scores)

    # pick most likely component
    labels = np.argmax(scores, axis=0)

    n = len(arr)
    out = []

    window_size = k + 1  # include itself

    # Slide window
    for i in range(n):
        target = arr[i]
        target_label = labels[i]

        l, r = i - 1, i + 1
        count = 1

        while count < window_size:
            left_ok = (l >= 0 and labels[l] == target_label)
            right_ok = (r < n and labels[r] == target_label)

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

        subset = np.array(arr[l + 1:r])

        if np.mean(subset) > 2.5 and np.var(subset) > 10:
            res = NN_ecDNA_predict(
                    model,
                    subset.tolist(),
                    out_log,
                    genes,
                    decision_rule,
                    verbose
                )[0]
            if verbose :
                out_log.append(f"{genes[0]} ecDNA:")
            out.append(res)
            
    return sum(out) > (ecDNA_percentage_thresh * n), sum(out)
