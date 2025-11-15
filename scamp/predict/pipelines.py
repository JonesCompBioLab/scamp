"""
A pipeline for predicting ecDNA status from single-cell copy-number
distributions.
"""

import numpy as np
import pandas as pd
import torch
import scanpy as sc

from scamp import io
from scamp import models
from scamp.predict import utilities
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster



def predict_ecdna_from_anndata(
    anndata_file,
    saved_model_directory,
    decision_rule,
    min_copy_number,
    max_percentile,
    filter_copy_number,
    cluster_distance_threshold
):
    counts_df = io.read_anndata_file(anndata_file)
    return predict(counts_df,
    saved_model_directory,
    decision_rule,
    min_copy_number,
    max_percentile,
    filter_copy_number,
    cluster_distance_threshold)


def predict_ecdna_from_mex(
    mex_folder,
    saved_model_directory,
    decision_rule,
    min_copy_number,
    max_percentile,
    filter_copy_number,
    cluster_distance_threshold
):
    counts_df = io.read_mex_file(mex_folder)

    return predict(counts_df,
    saved_model_directory,
    decision_rule,
    min_copy_number,
    max_percentile,
    filter_copy_number,
    cluster_distance_threshold)



def predict_ecdna_from_copy_number(
    counts_file,
    saved_model_directory,
    decision_rule,
    min_copy_number,
    max_percentile,
    filter_copy_number,
    cluster_distance_threshold
):

    counts_df = io.read_copy_numbers_file(counts_file)
    return predict(counts_df,
    saved_model_directory,
    decision_rule,
    min_copy_number,
    max_percentile,
    filter_copy_number,
    cluster_distance_threshold)


def predict(
    counts_df,
    saved_model_directory,
    decision_rule,
    min_copy_number,
    max_percentile,
    filter_copy_number,
    cluster_distance_threshold
) :
    model = models.SCAMP.load(saved_model_directory)

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

    prediction_df = cluster(prediction_df, counts_df, cluster_distance_threshold)

    return prediction_df

def cluster (
    prediction_df,
    counts_df,
    cluster_distance_threshold   
) :
    # Get each ecDNA's copy numbers as a vector
    ecDNA_genes = prediction_df.loc[prediction_df["pred"], "gene"].tolist()

    counts_df_ecDNA = counts_df[ecDNA_genes]
    gene_vectors = counts_df_ecDNA.T

    # Euclidean distance clustering
    Z = linkage(pdist(gene_vectors, metric="euclidean"), method="average")
    clusters = fcluster(Z, t=float(cluster_distance_threshold), criterion="distance")
    cluster_map = pd.Series(clusters, index=ecDNA_genes)

    # Add to dataframe
    prediction_df["cluster"] = -1
    prediction_df.loc[prediction_df["gene"].isin(ecDNA_genes), "cluster"] = clusters

    return prediction_df