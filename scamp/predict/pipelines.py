"""
A pipeline for predicting ecDNA status from single-cell copy-number
distributions.
"""

import numpy as np
import pandas as pd
import torch

from scamp import io
from scamp import models
from scamp.predict import utilities


def predict_ecdna_from_copy_number(
    counts_file,
    saved_model_directory,
    decision_rule,
    min_copy_number,
    max_percentile,
    filter_copy_number,
    whitelist_file
):

    counts_df = io.read_copy_numbers_file(counts_file)
    model = models.SCAMP.load(saved_model_directory)

    # subset by whitelist
    if whitelist_file:
        whitelist = pd.read_csv(whitelist_file, header=None).iloc[:,0].values
        counts_df = counts_df.loc[np.intersect1d(counts_df.index, whitelist)]

    X, genes_pass_filter = model.prepare_copy_numbers(
        counts_df.to_numpy(),
        np.array(counts_df.columns),
        min_copy_number=min_copy_number,
        max_percentile=max_percentile,
        filter_copy_number=filter_copy_number,
    )

    probas = model.proba(torch.Tensor(X)).detach().numpy()[:, 1]

    prediction_df = pd.DataFrame(X[:, 0:3])
    prediction_df.columns = ["mean", "var", "dispersion"]

    prediction_df["gene"] = genes_pass_filter
    prediction_df["proba"] = probas
    prediction_df["pred"] = prediction_df["proba"] >= decision_rule

    return prediction_df
