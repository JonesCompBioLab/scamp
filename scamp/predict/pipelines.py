"""
A pipeline for predicting ecDNA status from single-cell copy-number
distributions.
"""

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import os
import numpy as np
import pandas as pd
import torch
import scanpy as sc
import time
from pathlib import Path

from scamp import io
from scamp import models
from scamp.predict import utilities
from scamp import plotting



def predict_ecdna_from_anndata(
    out_log, anndata_file,
    saved_model_directory, whitelist_file,
    decision_rule,
    min_copy_number,
    max_percentile,
    filter_copy_number
):
    counts_df = io.read_anndata_file(anndata_file)
    return predict(out_log, counts_df,
    saved_model_directory, whitelist_file,
    decision_rule,
    min_copy_number,
    max_percentile,
    filter_copy_number)


def predict_ecdna_from_mex(
    out_log, mex_folder,
    saved_model_directory, whitelist_file,
    decision_rule,
    min_copy_number,
    max_percentile,
    filter_copy_number
):
    counts_df = io.read_mex_file(mex_folder)

    return predict(out_log, counts_df,
    saved_model_directory, whitelist_file,
    decision_rule,
    min_copy_number,
    max_percentile,
    filter_copy_number)



def predict_ecdna_from_copy_number(
    out_log, counts_file,
    saved_model_directory, whitelist_file,
    decision_rule,
    min_copy_number,
    max_percentile,
    filter_copy_number
):

    counts_df = io.read_copy_numbers_file(counts_file)
    return predict(out_log, counts_df,
    saved_model_directory, whitelist_file,
    decision_rule,
    min_copy_number,
    max_percentile,
    filter_copy_number)


def predict(
    out_log,
    counts_df,
    saved_model_directory,
    whitelist_file,
    decision_rule,
    min_copy_number,
    max_percentile,
    filter_copy_number
) :
    model = models.SCAMP.load(saved_model_directory)

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

def run_sample(file, output_dir, model_file, whitelist_file, decision_rule, min_copy_number, max_percentile, 
               filter_copy_number, no_plot) :
    out_log = []
    # Detect extension
    out_log.append(f'Running {file}')
    p = Path(file)
    if os.path.isdir(p) :
        has_matrix = any(p.glob("matrix.mtx*"))
        has_barcodes = any(p.glob("barcodes.tsv*"))
        if has_matrix and has_barcodes :
            mode = "MEX"
        else :
            out_log.append(f"{file} is non-MEX folder, skipping")
            return out_log

    else :
        copy_numbers_ext = file.split('.')[-1]
        if copy_numbers_ext == "h5ad" :
            mode = "anndata"
        elif copy_numbers_ext == 'csv' or copy_numbers_ext == 'tsv' :
            mode = "copynumber"
        else :
            out_log.append(f"{file} does not have extension h5ad, tsv, or csv. Skipping")
            return

    out_log.append(f"Running {file}")
    start = time.time()
    # Call different wrapper for each prediction type
    if mode == "copynumber":
        predictions = predict_ecdna_from_copy_number(
            out_log, file,
            model_file, whitelist_file,
            decision_rule,
            min_copy_number,
            max_percentile,
            filter_copy_number
        )
    elif mode == "MEX" :
        predictions = predict_ecdna_from_mex(
            out_log, file,
            model_file, whitelist_file,
            decision_rule,
            min_copy_number,
            max_percentile,
            filter_copy_number
        )
    else :
        predictions = predict_ecdna_from_anndata(
            out_log, file,
            model_file, whitelist_file,
            decision_rule,
            min_copy_number,
            max_percentile,
            filter_copy_number
        )

    os.makedirs(output_dir, exist_ok=True)

    # Output predictions and visualizations
    filename = Path(file).stem

    predictions.to_csv(f"{output_dir}/ecDNA_preds_{filename}.tsv", sep='\t')
    if not no_plot:
        plotting.plot_scamp_predictions_plotly(
            predictions,
            f"{output_dir}/ecDNA_predictions_{filename}.html",
            title=f"scAmp predictions for {filename.split('/')[-1]}"
        )

    end = time.time()


    out_log.append(f"File {file} completed in {end - start:.2f} seconds")
    return out_log

