"""
Functions to read in data.
"""

import numpy as np
import pandas as pd
import torch
import anndata as ad
from scipy import sparse



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
    ext = filename.split('.')[-1]
    if ext == "csv" :
        sep = ','
    elif ext == "tsv" :
        sep = '\t'

    counts_df = pd.read_csv(filename, sep=sep)

    return counts_df

def load_model(model_file):
    """Read in PyTorch model file.

    Args:
        model_file: File path to PyTorch model. 
    """

    model = torch.load(model_file, weights_only=False)
    return model

def read_anndata_file(filename, n_threads: int = None):
    """Read in copy number file.

    Reads in copy-number dataframe. Currently just is a wrapper for pd.DataFrame
    but we abstract this away in anticipation of doing other procedures
    on these copy-number files.

    Args:
        filename: Filename of copy-number data frame.
        n_threads: Number of threads to use. If None, number of physical cpu's
            of your system are used.
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