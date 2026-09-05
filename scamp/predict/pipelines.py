"""
A pipeline for predicting ecDNA status from single-cell copy-number
distributions.
"""

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import logging
import os
from functools import wraps
import inspect
from pprint import pformat
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


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def log_prediction_run(function):
    """Log the lifecycle of a prediction run to its output directory."""
    @wraps(function)
    def wrapper(file, output_dir, *args, **kwargs):
        os.makedirs(output_dir, exist_ok=True)
        file_handler = logging.FileHandler(Path(output_dir) / "scamp.log")
        file_handler.setFormatter(
            logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        )
        logger.addHandler(file_handler)
        start = time.perf_counter()

        try:
            parameters = inspect.signature(function).bind(
                file, output_dir, *args, **kwargs
            )
            parameters.apply_defaults()
            logger.info("Run parameters:\n%s", pformat(parameters.arguments))
            logger.info("Running prediction for %s", file)
            predictions = function(file, output_dir, *args, **kwargs)
        except Exception:
            logger.exception("Prediction failed for %s", file)
            raise
        else:
            logger.info("File %s completed in %.2f seconds", file, time.perf_counter() - start)
            return predictions
        finally:
            logger.removeHandler(file_handler)
            file_handler.close()

    return wrapper


def predict_ecdna_from_anndata(
    anndata_file: str,
    saved_model_directory: str,
    whitelist_file: str,
    decision_rule: float,
    min_copy_number: float,
    max_percentile: float,
    filter_copy_number: float
) -> pd.DataFrame:
    """Runs prediction pipeline from AnnData input.

    Args:
        anndata_file: Path to AnnData File
        saved_model_directory: Path to saved model directory
        whitelist_file: Path to whitelist file
        decision_rule: Likelihood decision rule
        min_copy_number: Minimum copy-number to consider
        max_percentile: Maximum percentile to cap copy-numbers
        filter_copy_number: Drop genes whose mean copy-number is below this threshold

    Returns:
        A pandas DataFrame with predictions and probabilities for each gene.
    """

    counts_df = io.read_anndata_file(anndata_file)
    predictions = _predict_ecDNA(counts_df,
        saved_model_directory, whitelist_file,
        decision_rule,
        min_copy_number,
        max_percentile,
        filter_copy_number)

    return predictions 


def predict_ecdna_from_mex(
    mex_folder: str,
    saved_model_directory: str,
    whitelist_file: str,
    decision_rule: float,
    min_copy_number: float,
    max_percentile: float,
    filter_copy_number: float
) -> pd.DataFrame:
    """Runs prediction pipeline from MEX input.

    Args:
        mex_folder: Path to MEX folder
        saved_model_directory: Path to saved model directory
        whitelist_file: Path to whitelist file
        decision_rule: Likelihood decision rule
        min_copy_number: Minimum copy-number to consider
        max_percentile: Maximum percentile to cap copy-numbers
        filter_copy_number: Drop genes whose mean copy-number is below this threshold

    Returns:
        A pandas DataFrame with predictions and probabilities for each gene.
    """

    counts_df = io.read_mex_file(mex_folder)

    predictions = _predict_ecDNA(counts_df,
        saved_model_directory, whitelist_file,
        decision_rule,
        min_copy_number,
        max_percentile,
        filter_copy_number)
    
    return predictions

def predict_ecdna_from_copy_number(
    counts_file: str,
    saved_model_directory: str,
    whitelist_file: str,
    decision_rule: float,
    min_copy_number: float,
    max_percentile: float,
    filter_copy_number: float
) -> pd.DataFrame:
    """Runs prediction pipeline from copy-number counts TSV.

    Args:
        counts_file: Path to copy-number counts file
        saved_model_directory: Path to saved model directory
        whitelist_file: Path to whitelist file
        decision_rule: Likelihood decision rule
        min_copy_number: Minimum copy-number to consider
        max_percentile: Maximum percentile to cap copy-numbers
        filter_copy_number: Drop genes whose mean copy-number is below this threshold

    Returns:
        A pandas DataFrame with predictions and probabilities for each gene.
    """

    counts_df = io.read_copy_numbers_file(counts_file)

    predictions = _predict_ecDNA(counts_df,
        saved_model_directory, whitelist_file,
        decision_rule,
        min_copy_number,
        max_percentile,
        filter_copy_number)
    
    return predictions

def _predict_ecDNA(
    counts_df: pd.DataFrame,
    saved_model_directory: str,
    whitelist_file: str,
    decision_rule: float,
    min_copy_number: float,
    max_percentile: float,
    filter_copy_number: float
) -> pd.DataFrame:
    """Runs prediction pipeline.

    Args:
        counts_df: A pandas DataFrame with copy-number counts
        saved_model_directory: Path to saved model directory
        whitelist_file: Path to whitelist file
        decision_rule: Likelihood decision rule
        min_copy_number: Minimum copy-number to consider
        max_percentile: Maximum percentile to cap copy-numbers
        filter_copy_number: Drop genes whose mean copy-number is below this threshold

    Returns:
        A pandas DataFrame with predictions and probabilities for each gene.
    """


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

@log_prediction_run
def predict_ecDNA_in_sample(
    file,
    output_dir,
    model_file,
    whitelist_file,
    decision_rule,
    min_copy_number,
    max_percentile, 
    filter_copy_number,
    no_plot,
) -> pd.DataFrame | None:
    """Runs the prediction pipeline on a single sample and outputs predictions and visualizations.

    Args:
        file: Path to input file (AnnData, MEX, or copy-number counts) 
        output_dir: Directory to save predictions and visualizations
        model_file: Path to saved model directory
        whitelist_file: Path to whitelist file
        decision_rule: Likelihood decision rule
        min_copy_number: Minimum copy-number to consider
        max_percentile: Maximum percentile to cap copy-numbers
        filter_copy_number: Drop genes whose mean copy-number is below this threshold
        no_plot: If True, do not generate visualizations
    
    Returns:
        A pandas DataFrame with predictions, or ``None`` when the input is skipped.
    """

    # Detect extension
    p = Path(file)
    if os.path.isdir(p):
        has_matrix = any(p.glob("matrix.mtx*"))
        has_barcodes = any(p.glob("barcodes.tsv*"))
        if has_matrix and has_barcodes:
            mode = "MEX"
        else:
            logger.warning("%s is a non-MEX folder; skipping", file)
            return None

    else:
        copy_numbers_ext = file.split('.')[-1]
        if copy_numbers_ext == "h5ad":
            mode = "anndata"
        elif copy_numbers_ext == "csv" or copy_numbers_ext == "tsv":
            mode = "copynumber"
        else:
            logger.warning(
                "%s does not have extension h5ad, tsv, or csv; skipping", file
            )
            return None

    # Call different wrapper for each prediction type
    if mode == "copynumber":
        predictions = predict_ecdna_from_copy_number(
            file,
            model_file, whitelist_file,
            decision_rule,
            min_copy_number,
            max_percentile,
            filter_copy_number
        )
    elif mode == "MEX":
        predictions = predict_ecdna_from_mex(
            file,
            model_file, whitelist_file,
            decision_rule,
            min_copy_number,
            max_percentile,
            filter_copy_number
        )
    else :
        predictions = predict_ecdna_from_anndata(
            file,
            model_file, whitelist_file,
            decision_rule,
            min_copy_number,
            max_percentile,
            filter_copy_number
        )

    os.makedirs(output_dir, exist_ok=True)

    # Output predictions and visualizations
    filename = Path(file).stem

    predictions.to_csv(f"{output_dir}/ecDNA_predictions_{filename}.tsv", sep='\t')

    if not no_plot:
        plotting.plot_scamp_predictions_plotly(
            predictions,
            f"{output_dir}/ecDNA_predictions_{filename}.html",
            title=f"scAmp predictions for {filename.split('/')[-1]}"
        )

    return predictions

