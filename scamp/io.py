"""
Functions to read in data.
"""

import numpy as np
import pandas as pd
import torch


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

    counts_df = pd.read_csv(filename, sep='\t', index_col=0)

    return counts_df

def load_model(model_file):
    """Read in PyTorch model file.

    Args:
        model_file: File path to PyTorch model. 
    """

    model = torch.load(model_file, weights_only=False)
    return model