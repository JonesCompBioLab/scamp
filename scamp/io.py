"""
Functions to read in data.
"""

import os
from pathlib import Path
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import numpy as np
import pandas as pd
import torch
import anndata as ad
from scipy import sparse
import scanpy as sc

from scamp.mixins.errors import CLIError


def detect_copy_number_file_type(filename):
    """Detect copy number file type.

    Args:
        filename: Filename of copy-number data frame.
    """

    p = Path(filename)
    mode = None
    if os.path.isdir(p):
        has_matrix = any(p.glob("matrix.mtx*"))
        has_barcodes = any(p.glob("barcodes.tsv*"))
        if has_matrix and has_barcodes:
            mode = "MEX"
        else:
            # out_log.append(f"{file} is non-MEX folder, skipping")
            raise CLIError("Detected non-MEX folder, exiting.")

    else:
        copy_numbers_ext = filename.split(".")[-1]
        if copy_numbers_ext == "h5ad":
            mode = "anndata"
        elif copy_numbers_ext == "csv" or copy_numbers_ext == "tsv":
            mode = "text"
        else:
            # out_log.append(f"{file} does not have extension h5ad, tsv, or csv. Skipping")
            raise CLIError("Detected non-supported file type: does not have "
                           "extension h5ad, tsv or csv, exiting.")

    return mode

def read_copy_numbers_file(filename, n_threads: int = None):
    """Read in copy number file.

    Reads in copy-number dataframe. Currently just is a wrapper for pd.DataFrame
    but we abstract this away in anticipation of doing other procedures
    on these copy-number files.

    Args:
        filename: Filename of copy-number data frame.
        n_threads: Number of threads to use. If None, number of physical cpu's
            of your system are used.
    """

    mode = detect_copy_number_file_type(filename)

    if mode == "text":
        ext = filename.split('.')[-1]
        if ext == "csv" :
            sep = ','
        else :
            sep = '\t'

        counts_df = pd.read_csv(filename, sep=sep, index_col=0)

    elif mode == 'anndata':
        counts_df = read_anndata_file(filename, n_threads)

    elif mode == 'MEX':
        counts_df = read_mex_file(filename, n_threads)

    return counts_df


def read_anndata_file(filename):
    """Read in anndata file.

    Reads in h5ad

    Args:
        filename: Filename of anndata h5ad.

    Returns:
        counts_df: DataFrame of counts with cells as rows and genes as columns.

    """
    adata = ad.read_h5ad(filename)
    X = adata.X
    if sparse.issparse(X):
        X = X.toarray()

    counts_df = pd.DataFrame(
        X,
        index=adata.obs_names,
        columns=adata.var_names
    )

    return counts_df

def read_mex_file(filename):
    """Read in anndata file.

    Reads in h5ad

    Args:
        filename: Filename of anndata h5ad.
    
    Returns:
        counts_df: DataFrame of counts with cells as rows and genes as columns.
    """
    adata = sc.read_10x_mtx(filename)

    X = adata.X
    if sparse.issparse(X):
        X = X.toarray()

    counts_df = pd.DataFrame(
        X,
        index=adata.obs_names,
        columns=adata.var_names
    )
    return counts_df


def load_model(model_file):
    """Read in PyTorch model file.

    Args:
        model_file: File path to PyTorch model. 
    """

    model = torch.load(model_file, weights_only=False)
    return model