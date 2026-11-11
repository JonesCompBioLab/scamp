import heapq
import os
import shutil
from collections import defaultdict
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from cnmf import cNMF
from cnmf.cnmf import efficient_ols_all_cols
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.optimize import nnls
from sklearn.metrics import pairwise_distances, silhouette_score


def remove_duplicates(cellbygene, tol=1e-5):
    """Collapses genes that have the same copy-number profile across cells.

    Args:
        cellbygene: Cell-by-gene copy-number DataFrame.
        tol: Largest absolute difference, in every cell, for two genes to count
            as duplicates.

    Returns:
        A tuple `(df_unique, column_groups)`. `df_unique` keeps the first gene of
        each group of duplicates. `column_groups` maps each kept gene to all the
        genes it stands for, including itself.
    """
    keep_cols = []
    column_groups = {}

    for col in cellbygene.columns:
        found_match = False

        for rep in keep_cols:
            if np.allclose(cellbygene[col], cellbygene[rep], atol=tol, rtol=0):
                column_groups[rep].append(col)
                found_match = True
                break

        if not found_match:
            keep_cols.append(col)
            column_groups[col] = [col]

    df_unique = cellbygene[keep_cols]
    return df_unique, column_groups


def cNMF_deconvolution(
    cellbygene_df,
    sample_name,
    out_log,
    n_iter=50,
    error_w=0.25,
    score_cutoff=3,
    ecDNA_remove_cutoff=0.01,
    log_dir="./temp",
    density_threshold=0.1,
    hier_ddist=1.3,
    max_species=None,
    num_ecDNA=None,
    seed=None,
    verbose=False,
):
    """Groups genes into ecDNA species with cNMF, assigning genes by ablation.

    cNMF is run on copy-numbers above the chromosomal baseline of 2. Each gene's
    weight in each species is then set to zero in turn, and the resulting increase
    in reconstruction error is that gene's score for the species. A gene joins
    every species where its score is above `score_cutoff`, or the species with its
    highest score if there are none. Species left with no genes are dropped.

    Writes intermediate files to `log_dir` and removes the cNMF temporary folder
    when done.

    Args:
        cellbygene_df: Cell-by-gene copy-number DataFrame, still including the
            chromosomal baseline of 2.
        sample_name: Name for the cNMF run and its log folder.
        out_log: List of log lines; cNMF's output is appended to it.
        n_iter: Number of NMF runs per number of species.
        error_w: Weight of reconstruction error when choosing the number of species.
        score_cutoff: Score above which a gene is assigned to a species.
        ecDNA_remove_cutoff: Species whose highest spectra score is below this are
            dropped, unless that would drop all of them.
        log_dir: Directory for cNMF intermediate files.
        density_threshold: Density threshold for the cNMF consensus.
        hier_ddist: Distance for null nodes in hierarchical clustering, used to set
            `max_species` when it is None.
        max_species: Largest number of species to check. If None, it is the number
            chosen by hierarchical clustering.
        num_ecDNA: Number of species, if known. If None, it is chosen from the
            cNMF stability and reconstruction error.
        seed: Random seed for cNMF.
        verbose: Whether to log details of the deconvolution.

    Returns:
        A tuple `(species_to_gene, cell_by_ecDNA)`. `species_to_gene` maps each
        species name ("pred_ecDNA_<n>", numbered from 1 after dropped species are
        removed) to its genes, and `cell_by_ecDNA` is a cell-by-species DataFrame of
        ecDNA amounts estimated with NNLS, with the same species as columns.
    """
    if num_ecDNA is not None and max_species is not None and num_ecDNA > max_species:
        raise ValueError(
            f"num_ecDNA ({num_ecDNA}) is larger than max_species ({max_species})."
        )

    os.makedirs(log_dir, exist_ok=True)

    # Decrease all by 2, and make none less than zero
    cellbygene = cellbygene_df - 2
    cellbygene = cellbygene.clip(lower=0)

    # cNMF raises "N cells have zero counts of overdispersed genes ... Quitting!" if any
    # cell is all-zero here, i.e. has no copy-number above 2 at any ecDNA+ gene. Such
    # cells are left out of the cNMF input and get 0 usage afterwards, so the output
    # keeps all cells (ported from the old local fix U7).
    nonzero = cellbygene.sum(axis=1) > 0
    n_dropped = int((~nonzero).sum())
    if n_dropped > 0:
        out_log.append(
            f"cNMF: left out {n_dropped} of {len(cellbygene)} cells with no copy-number "
            f"above 2 at any ecDNA+ gene; added back with 0 usage"
        )

    # Export (needs to be an outputted csv otherwise we cannot)
    cellbygene_path = f"{log_dir}/cellbygene_minus_chromosomal.tsv"
    cellbygene[nonzero].to_csv(cellbygene_path, sep="\t")

    if num_ecDNA is None:
        num_ecDNA = _cNMF_find_k(
            max_species,
            cellbygene_df,
            hier_ddist,
            log_dir,
            sample_name,
            cellbygene_path,
            n_iter,
            seed,
            error_w,
            out_log,
            verbose,
        )

    # Silence certain printouts
    out = StringIO()
    err = StringIO()

    with redirect_stdout(out), redirect_stderr(err):

        cnmf_obj = cNMF(output_dir=log_dir, name=sample_name)
        cnmf_obj.prepare(
            counts_fn=cellbygene_path,
            tpm_fn=cellbygene_path,
            components=num_ecDNA,
            n_iter=n_iter,
            seed=seed,
        )

        input_counts = pd.read_csv(cellbygene_path, sep="\t", index_col=0)
        adata = ad.AnnData(input_counts)
        cnmf_obj.save_norm_counts(adata)
        cnmf_obj.factorize(worker_i=0, total_workers=1)
        cnmf_obj.combine()

        cnmf_obj.consensus(
            k=num_ecDNA,
            density_threshold=density_threshold,
            close_clustergram_fig=True,
            refit_usage=False,
        )
        usage_df, spectra_scores, spectra_tpm, top_genes = cnmf_obj.load_results(
            K=num_ecDNA, density_threshold=density_threshold, norm_usage=False
        )

        rf_usages = pd.read_csv(
            f"{log_dir}/{sample_name}/{sample_name}.usages.k_{num_ecDNA}.dt_{str(density_threshold).replace('.', '_')}.consensus.txt",
            sep="\t",
            index_col=0,
        )

        # Cells left out of cNMF have 0 usage; everything below uses all cells
        assert list(rf_usages.index) == list(cellbygene.index[nonzero])
        rf_usages = rf_usages.reindex(cellbygene.index, fill_value=0.0)
        usage_df = usage_df.reindex(cellbygene.index, fill_value=0.0)

        Y = cellbygene.values
        unnormalized_spectra_scores = efficient_ols_all_cols(rf_usages.values, Y, normalize_y=False)

    out_log.append(out.getvalue())
    out_log.append(err.getvalue())

    # [ecDNA, gene, iter]
    rows, cols = unnormalized_spectra_scores.shape
    beta2 = unnormalized_spectra_scores.copy()

    base_loss = np.linalg.norm(Y - (rf_usages.values @ unnormalized_spectra_scores), "fro") ** 2 / (
        Y.shape[0] * Y.shape[1]
    )
    if verbose:
        out_log.append(f"Base loss: {base_loss}")

    losses = np.zeros((rows, cols))

    for i in range(rows):
        for j in range(cols):
            temp = beta2[i, j]
            beta2[i, j] = 0
            usages_new = cellbygene.values @ np.linalg.pinv(beta2)
            Y_hat = usages_new @ beta2
            loss = np.linalg.norm(Y - Y_hat, "fro") ** 2 / (Y.shape[0] * Y.shape[1])
            beta2[i, j] = temp
            losses[i, j] = loss - base_loss

    unnormalized_spectra_scores = pd.DataFrame(
        unnormalized_spectra_scores.T, columns=rf_usages.columns, index=cellbygene.columns
    )
    losses = pd.DataFrame(losses.T, columns=rf_usages.columns, index=cellbygene.columns)

    # Parse cNMF spectra scores
    unnormalized_spectra_scores.columns = [
        f"pred_ecDNA_{col}" for col in unnormalized_spectra_scores.columns
    ]
    losses.columns = [f"pred_ecDNA_{col}" for col in losses.columns]
    spectra_scores.columns = [f"pred_ecDNA_{col}" for col in spectra_scores.columns]
    spectra_tpm.columns = [f"pred_ecDNA_{col}" for col in spectra_tpm.columns]

    # Remove ecDNA that don't pass threshold
    cols_to_remove = spectra_scores.columns[spectra_scores.max() < ecDNA_remove_cutoff].tolist()
    if len(cols_to_remove) < len(spectra_scores.columns):
        unnormalized_spectra_scores = unnormalized_spectra_scores.drop(columns=cols_to_remove)

    ecDNA_species = unnormalized_spectra_scores.columns
    species_to_gene = {}

    # Keep track of maximum value for the gene
    used_genes = set()
    all_genes = set()
    max_species_gene = {}
    max_species_gene_value = {}

    for species in ecDNA_species:
        species_to_gene[species] = []
    for i, row in losses.iterrows():
        all_genes.add(i)
        for species in ecDNA_species:
            if row[species] > score_cutoff:
                species_to_gene[species].append(i)
                used_genes.add(i)
            if i in max_species_gene_value:
                if row[species] > max_species_gene_value[i]:
                    max_species_gene_value[i] = row[species]
                    max_species_gene[i] = species
            else:
                max_species_gene_value[i] = row[species]
                max_species_gene[i] = species

    # Put all the ones who aren't in any one in an ecDNA
    for empty_gene in all_genes - used_genes:
        if empty_gene in max_species_gene:
            species_to_gene[max_species_gene[empty_gene]].append(empty_gene)
        else:
            if verbose:
                out_log.append(f"Warning: {empty_gene} not in any species")

    # Remove components with no genes (likely batch effects)
    for species in ecDNA_species:
        if len(species_to_gene[species]) == 0:
            del species_to_gene[species]
    ecDNA_species = species_to_gene.keys()

    cell_by_ecDNA = _cNMF_usage(ecDNA_species, species_to_gene, spectra_tpm, cellbygene, usage_df)

    # Number the kept species 1, 2, ... in cNMF's order, so the names don't skip the
    # programs removed above, and keep only their usage columns
    new_names = {old: f"pred_ecDNA_{i + 1}" for i, old in enumerate(species_to_gene)}
    if verbose:
        out_log.append(f"cNMF species renamed: {new_names}")
    species_to_gene = {new_names[old]: genes for old, genes in species_to_gene.items()}
    cell_by_ecDNA = cell_by_ecDNA[list(new_names)].rename(columns=new_names)

    shutil.rmtree(f"{log_dir}/{sample_name}/cnmf_tmp")

    return species_to_gene, cell_by_ecDNA


def residualize(Y, X):
    """Removes the linear effect of `X` from `Y`.

    Args:
        Y: Array of shape (n_cells, n_targets) or (n_cells,).
        X: Covariates, of shape (n_cells, n_covariates).

    Returns:
        The residuals of a least-squares fit of `Y` on `X`, with the same shape as `Y`.
    """
    XtX = X.T @ X
    XtX_inv = np.linalg.inv(XtX + 1e-12 * np.eye(X.shape[1]))
    P = X @ XtX_inv @ X.T
    return Y - P @ Y


def fwl_beta(Xk, Y_res):
    """Least-squares coefficient of `Y_res` on a single covariate `Xk`.

    With `Xk` and `Y_res` both residualized on the other covariates (see
    `residualize`), this is that covariate's coefficient in the full regression,
    by the Frisch-Waugh-Lovell theorem.

    Args:
        Xk: Covariate, of shape (n_cells,).
        Y_res: Residualized targets, of shape (n_cells, n_targets) or (n_cells,).

    Returns:
        The coefficient for each target, of shape (n_targets,), or a scalar.
    """
    return (Xk.T @ Y_res) / (Xk.T @ Xk + 1e-12)


def hierarchical_deconvolution(cellbygene_df, hier_ddist=1.3, max_species=None, num_ecDNA=None):
    """Groups genes into ecDNA species by hierarchical clustering of their correlations.

    Args:
        cellbygene_df: Cell-by-gene copy-number DataFrame, still including the
            chromosomal baseline of 2.
        ddist: Distance for the null node in the silhouette score used to choose
            the number of species.
        max_species: Largest number of species to check. If None, checks up to the
            number of genes minus 1.
        num_ecDNA: Number of species, if known. If None, it is chosen by silhouette score.

    Returns:
        A tuple `(species_to_gene, cell_by_ecDNA)`. `species_to_gene` maps each
        species name ("pred_ecDNA_<n>") to its genes. `cell_by_ecDNA` is a
        cell-by-species DataFrame: for each species, the mean copy-number above 2
        of its genes whose total is within 1.3 times that of its lowest gene.
    """

    X = cellbygene_df
    embed = np.corrcoef(X, rowvar=False)

    if num_ecDNA is None:
        num_ecDNA = _hier_get_k(embed, max_species, hier_ddist)

    Z = linkage(embed, method="average", metric="correlation")
    clusters = fcluster(Z, t=num_ecDNA, criterion="maxclust")

    species_to_gene = defaultdict(list)
    for i in range(len(clusters)):
        species_to_gene[f"pred_ecDNA_{clusters[i]}"].append(cellbygene_df.columns[i])
    reversed_observed = defaultdict(list)

    for key, values in species_to_gene.items():
        for v in values:
            reversed_observed[v].append(key)

    # When calculating usage do subtract 2
    cellbygene_temp = cellbygene_df - 2
    cellbygene_temp = cellbygene_temp.clip(lower=0)

    cellbyecDNA = pd.DataFrame(0.0, index=cellbygene_temp.index, columns=species_to_gene.keys())

    for species in list(species_to_gene.keys()):
        # Just trust the extra counts of the smallest one and those 1.3 times at most above it (which does not have duplicates hopefully or is on multiple ecDNA)
        genes = species_to_gene[species]
        gene_sums = cellbygene_temp[genes].sum()
        min_value = gene_sums.min()
        threshold = 1.3 * min_value
        genes_within_range = gene_sums[gene_sums <= threshold].index.tolist()
        subset = cellbygene_temp[genes_within_range]
        avg_list = subset.mean(axis=1).values
        cellbyecDNA[species] = avg_list

    return species_to_gene, cellbyecDNA


def hybrid_deconvolution(
    cellbygene_df,
    out_log,
    cNMF_thresh=0.55,
    sample_name="cNMF_sample",
    max_species=None,
    n_iter=50,
    error_w=0.25,
    score_cutoff=3,
    ecDNA_remove_cutoff=0.01,
    log_dir="./temp",
    density_threshold=0.1,
    hier_ddist=1.3,
    num_ecDNA=None,
    seed=None,
    verbose=False,
):
    """Groups genes into ecDNA species, choosing between hierarchical clustering and cNMF.

    Genes are first clustered hierarchically. If the clusters look like they
    overlap (genes shared between species), cNMF is run instead, with the
    hierarchical number of species as its upper bound. Otherwise the hierarchical
    clusters are kept, and their per-cell amounts are estimated with NNLS.

    Args:
        cellbygene_df: Cell-by-gene copy-number DataFrame, still including the
            chromosomal baseline of 2.
        out_log: List of log lines; messages are appended to it.
        cNMF_thresh: Each cluster's mean profile is fit as a combination of the
            others. Overlap is detected, and cNMF used, when the three lowest fit
            errors are all below this. Needs at least 3 clusters.
        sample_name: Name for the cNMF run and its log folder.
        max_species: Largest number of species to check when clustering
            hierarchically. If None, checks up to the number of genes minus 1.
        n_iter: Number of NMF runs per number of species.
        error_w: Weight of reconstruction error when choosing the number of cNMF species.
        score_cutoff: Score above which a gene is assigned to a cNMF species.
        ecDNA_remove_cutoff: cNMF species whose highest spectra score is below this
            are dropped.
        log_dir: Directory for cNMF intermediate files.
        density_threshold: Density threshold for the cNMF consensus.
        hier_ddist: Distance for null nodes in hierarchical clustering.
        num_ecDNA: Number of species, if known, used by both hierarchical clustering
            and cNMF. If None, hierarchical clustering chooses it by silhouette score,
            and cNMF chooses its own number up to that one.
        seed: Random seed for cNMF.
        verbose: Whether to log details of the deconvolution.

    Returns:
        A tuple `(species_to_gene, cell_by_ecDNA)`. `species_to_gene` maps each
        species name ("pred_ecDNA_<n>") to its genes, and `cell_by_ecDNA` is a
        cell-by-species DataFrame of estimated ecDNA amounts.
    """

    X = cellbygene_df
    embed = np.corrcoef(X, rowvar=False)

    # Number of hierarchical species; kept apart from num_ecDNA so cNMF still chooses
    # its own number (up to this one) when num_ecDNA is not given
    if num_ecDNA is None:
        hier_num_ecDNA = _hier_get_k(embed, max_species, hier_ddist)
    else:
        hier_num_ecDNA = num_ecDNA

    Z = linkage(embed, method="average", metric="correlation")
    clusters = fcluster(Z, t=hier_num_ecDNA, criterion="maxclust")

    # Continue with hierarchical
    if not _check_overlap(X, clusters, cNMF_thresh):
        return _hier_nnls(
            cellbygene_df, ddist=hier_ddist, max_species=max_species, num_ecDNA=hier_num_ecDNA
        )

    else:
        if verbose:
            out_log.append("Detected Overlap")
        return cNMF_deconvolution(
            cellbygene_df,
            sample_name,
            out_log,
            max_species=hier_num_ecDNA,
            n_iter=n_iter,
            error_w=error_w,
            score_cutoff=score_cutoff,
            ecDNA_remove_cutoff=ecDNA_remove_cutoff,
            log_dir=log_dir,
            density_threshold=density_threshold,
            hier_ddist=hier_ddist,
            num_ecDNA=num_ecDNA,
            seed=seed,
            verbose=verbose,
        )


def naive_deconvolution(cellbygene_df, threshold=0.65, num_ecDNA=None):
    """Groups genes into ecDNA species by hierarchical clustering of their copy-number profiles.

    Args:
        cellbygene_df: Cell-by-gene copy-number DataFrame, still including the
            chromosomal baseline of 2.
        threshold: Correlation distance at which to cut the clustering tree, used
            when `num_ecDNA` is None.
        num_ecDNA: Number of species, if known.

    Returns:
        A tuple `(species_to_gene, cell_by_ecDNA)`. `species_to_gene` maps each
        species name ("pred_ecDNA_<n>") to its genes. `cell_by_ecDNA` is a
        cell-by-species DataFrame: for each species, the mean copy-number above 2
        of its genes whose total is within 1.3 times that of its lowest gene.
    """
    X = cellbygene_df.T

    if num_ecDNA is None:
        Z = linkage(X, method="average", metric="correlation")
        clusters = fcluster(Z, t=threshold, criterion="distance")
        num_ecDNA = np.unique(clusters).size

    else:
        Z = linkage(X, method="average", metric="correlation")
        clusters = fcluster(Z, t=num_ecDNA, criterion="maxclust")

    species_to_gene = defaultdict(list)
    for i in range(len(clusters)):
        species_to_gene[f"pred_ecDNA_{clusters[i]}"].append(cellbygene_df.columns[i])
    reversed_observed = defaultdict(list)

    for key, values in species_to_gene.items():
        for v in values:
            reversed_observed[v].append(key)

    # When calculating usage do subtract 2
    cellbygene_temp = cellbygene_df - 2
    cellbygene_temp = cellbygene_temp.clip(lower=0)

    cellbyecDNA = pd.DataFrame(0.0, index=cellbygene_temp.index, columns=species_to_gene.keys())

    for species in list(species_to_gene.keys()):
        # Just trust the extra counts of the smallest one and those 1.3 times at most above it (which does not have duplicates hopefully or is on multiple ecDNA)
        genes = species_to_gene[species]
        gene_sums = cellbygene_temp[genes].sum()
        min_value = gene_sums.min()
        threshold = 1.3 * min_value
        genes_within_range = gene_sums[gene_sums <= threshold].index.tolist()
        subset = cellbygene_temp[genes_within_range]
        avg_list = subset.mean(axis=1).values
        cellbyecDNA[species] = avg_list

    return species_to_gene, cellbyecDNA


def _cNMF_find_k(
    max_species,
    cellbygene_df,
    hier_ddist,
    log_dir,
    sample_name,
    cellbygene_path,
    n_iter,
    seed,
    error_w,
    out_log,
    verbose,
):
    if max_species is None:
        embed = np.corrcoef(cellbygene_df, rowvar=False)
        max_species = _hier_get_k(embed, None, hier_ddist)

    # Silence certain printouts
    out = StringIO()
    err = StringIO()
    with redirect_stdout(out), redirect_stderr(err):

        cnmf_obj = cNMF(output_dir=log_dir, name=sample_name)
        check_one = False
        # Goes to one above because the next layer is needed for the score comparison
        counts_to_check = range(1, max_species + 2)
        if 1 in counts_to_check:
            cnmf_obj.prepare(
                counts_fn=cellbygene_path,
                tpm_fn=cellbygene_path,
                components=counts_to_check[1:],
                n_iter=n_iter,
                seed=seed,
            )
            check_one = True

        else:
            cnmf_obj.prepare(
                counts_fn=cellbygene_path,
                tpm_fn=cellbygene_path,
                components=counts_to_check,
                n_iter=n_iter,
                seed=seed,
            )

        input_counts = pd.read_csv(cellbygene_path, sep="\t", index_col=0)
        adata = ad.AnnData(input_counts)
        cnmf_obj.save_norm_counts(adata)
        cnmf_obj.factorize(worker_i=0, total_workers=1)
        cnmf_obj.combine()
        cnmf_obj.k_selection_plot(close_fig=True)

        # Find best number of ecDNA using stability and error

        npz = np.load(
            f"{log_dir}/{sample_name}/{sample_name}.k_selection_stats.df.npz", allow_pickle=True
        )

        k_df = pd.DataFrame(data=npz["data"], index=npz["index"], columns=npz["columns"])

        # Include 1 (stability always at 1)
        if check_one:
            cnmf_obj_1 = cNMF(output_dir=log_dir, name=sample_name)

            cnmf_obj_1.prepare(
                counts_fn=cellbygene_path,
                tpm_fn=cellbygene_path,
                components=1,
                n_iter=n_iter,
                seed=seed,
            )
            cnmf_obj_1.save_norm_counts(adata)
            with redirect_stdout(out), redirect_stderr(err):
                cnmf_obj_1.factorize(worker_i=0, total_workers=1)
            cnmf_obj_1.combine()

            norm_counts = sc.read(cnmf_obj_1.paths["normalized_counts"])

            with np.load(cnmf_obj_1.paths["merged_spectra"] % 1, allow_pickle=True) as f:
                obj = pd.DataFrame(**f)
                spectra = obj
            l2_spectra = (spectra.T / np.sqrt((spectra**2).sum(axis=1))).T
            median_spectra = pd.DataFrame(l2_spectra.median(axis=0)).T
            median_spectra = (median_spectra.T / median_spectra.sum(1)).T
            rf_usages = cnmf_obj_1.refit_usage(norm_counts.X, median_spectra)
            rf_usages = pd.DataFrame(rf_usages, index=norm_counts.obs.index)
            rf_pred = rf_usages.dot(median_spectra)

            if sp.issparse(norm_counts.X):
                prediction_error = ((norm_counts.X.todense() - rf_pred) ** 2).sum().sum()
            else:
                prediction_error = ((norm_counts.X - rf_pred) ** 2).sum().sum()

            new_row = pd.DataFrame(
                [
                    {
                        "k": 1,
                        "local_density_threshold": 0.5,
                        "silhouette": 1,
                        "prediction_error": prediction_error,
                    }
                ]
            )
            k_df = pd.concat([new_row, k_df], ignore_index=True)

    max_score = max(k_df["prediction_error"])

    k_df["normalized_prediction_error"] = (k_df["prediction_error"]) / (max_score)

    # Ratio between current prediction and the prediction score of k + 1 (minus 1, otherwise they are all greater than 1 really as it should decrease)
    # Default to something that can never be picked (realistically)
    k_df["prediction_ratio"] = 100 / error_w

    for count in counts_to_check:
        if count <= max_species:
            row_k = k_df[k_df["k"] == count].iloc[0]
            row_kp1 = k_df[k_df["k"] == count + 1].iloc[0]
            k_df.loc[k_df["k"] == count, "prediction_ratio"] = (
                (row_k["normalized_prediction_error"] + 1e-5)
                / (row_kp1["normalized_prediction_error"] + 1e-5)
            ) - 1

    k_df["score"] = k_df["silhouette"] - error_w * k_df["prediction_ratio"]
    num_ecDNA = int(k_df.loc[k_df["score"].idxmax()]["k"])

    if verbose:
        out_log.append(k_df.to_string())
        out_log.append(f"Number of ecDNA chosen: {num_ecDNA}")
    return num_ecDNA


# Version of hierarchical that uses nnls (best for hybrid since theoretically we know there are uniques)
def _hier_nnls(cellbygene_df, ddist=1.3, max_species=None, num_ecDNA=None):

    X = cellbygene_df
    embed = np.corrcoef(X, rowvar=False)

    if num_ecDNA is None:
        num_ecDNA = _hier_get_k(embed, max_species, ddist)

    Z = linkage(embed, method="average", metric="correlation")
    clusters = fcluster(Z, t=num_ecDNA, criterion="maxclust")

    species_to_gene = defaultdict(list)
    for i in range(len(clusters)):
        species_to_gene[f"pred_ecDNA_{clusters[i]}"].append(cellbygene_df.columns[i])
    reversed_observed = defaultdict(list)

    for key, values in species_to_gene.items():
        for v in values:
            reversed_observed[v].append(key)

    # When calculating usage do subtract 2
    cellbygene_temp = cellbygene_df - 2
    cellbygene_temp = cellbygene_temp.clip(lower=0)

    cellbyecDNA = pd.DataFrame(0.0, index=cellbygene_temp.index, columns=species_to_gene.keys())

    # NNLS setup: create species by gene matrix (assume no overlaps at this point)
    species_profiles = {}

    for pred_species, genes in species_to_gene.items():

        gene_sums = cellbygene_temp[genes].sum()

        min_value = gene_sums.min()
        threshold = 1.2 * min_value

        # These give us the average to divide by to get genes per species
        baseline_genes = gene_sums[gene_sums <= threshold]
        baseline_mean = baseline_genes.mean()

        profile = pd.Series(0.0, index=cellbygene_temp.columns, dtype=float)

        for gene in genes:

            profile[gene] = gene_sums[gene] / baseline_mean

        species_profiles[pred_species] = profile

    # gene by species matrix
    gene_by_species = pd.DataFrame(species_profiles)
    A = gene_by_species.values

    # Run NNLS
    pred_species_usage = []
    for cell in cellbygene_temp.index:
        b = cellbygene_temp.loc[cell].values
        x, residual = nnls(A, b)
        pred_species_usage.append(x)

    cellbyecDNA = pd.DataFrame(
        pred_species_usage, index=cellbygene_temp.index, columns=gene_by_species.columns
    )

    return species_to_gene, cellbyecDNA


# Gets the usage matrix for cNMF
def _cNMF_usage(ecDNA_species, species_to_gene, spectra_tpm, cellbygene, usage_df):
    # Find how to rescale usage in terms of tpm
    # Assumes the lowest of the used genes is 1, and takes the average of the lowest and those 1.2 times away (incase there are slight deviations)
    usage_scale = {}
    for species in ecDNA_species:
        usage_scale[species] = 1
        obs_genes = species_to_gene[species]
        if len(obs_genes) == 0:
            continue

        tpm_values = spectra_tpm.loc[obs_genes, species].values
        min_tpm = np.min(tpm_values)
        near_min = tpm_values[tpm_values <= min_tpm * 1.2]
        usage_scale[species] = np.mean(near_min)

    # Now recreate a better spectra_tpm
    for species in ecDNA_species:
        scale = usage_scale[species]

        if scale == 0:
            continue

        # divide by the found scale
        spectra_tpm[species] /= scale

        # zero out genes not in the ecDNA species
        mask = ~spectra_tpm.index.isin(species_to_gene[species])
        spectra_tpm.loc[mask, species] = 0

    # Run NNLS optimization to find better usage matrix
    def NNLS(cellbygene_df, spectra_tpm):
        X = cellbygene_df.values
        H = spectra_tpm.values.T

        W = np.zeros((X.shape[0], H.shape[0]))

        for i in range(X.shape[0]):
            W[i], _ = nnls(H.T, X[i])
        return W

    usage_df_new_vals = NNLS(cellbygene, spectra_tpm)
    usage_df.loc[:, :] = usage_df_new_vals
    # NNLS gives one usage column per spectra_tpm column, so they share its names
    usage_df.columns = spectra_tpm.columns

    return usage_df


# Silhouette score using an extra node with distance ddist to all others
def _silhouette_with_extra(X, labels, ddist, metric="correlation", eps=1e-5):
    n = X.shape[0]

    # Square the distances to slighly reduce the effect of double-correlation
    X = np.sign(X) * (X**2)

    # Add new dummy node to help with the correlations if everything is correlated
    X_new = np.zeros((n + 1, n + 1))
    X_new[:n, :n] = X
    X_new[n, n] = 1

    D = pairwise_distances(X_new, metric=metric)

    # Reset the distance so it doesn't overblow the results for k = 1
    D += eps
    D[n, :] = ddist
    D[:, n] = ddist

    np.fill_diagonal(D, 0)

    # Add extra label for the dummy node
    labels = np.append(labels, max(labels) + 1)

    # pass it as precomputed
    silhouette = silhouette_score(D, labels, metric="precomputed")

    return silhouette


# Calculate the best number of ecDNA species to try for hierarchical clustering
# Uses hierarchical clustering
# Returns the number of ecDNA
def _hier_get_k(X, max_species, ddist, leeway=0):
    best_num = -1
    best_silhouette = -1
    silhouettes = []

    if max_species is None:
        max_species, _ = X.shape
        max_species = max_species - 1

    nums_to_try = [i for i in range(1, max_species + 1)]

    for i in nums_to_try:
        Z = linkage(X, method="average", metric="euclidean")
        if i == 1:
            clusters = np.full((X.shape[0]), 1)
        else:
            clusters = fcluster(Z, t=i, criterion="maxclust")

        silhouette = _silhouette_with_extra(X, clusters, ddist, metric="euclidean")
        silhouettes.append(silhouette)
        if silhouette > best_silhouette:
            best_silhouette = silhouette

    # Allow some leeway around the silhouette score, to favor greater values wiht just slightly worse silhouette scores
    for idx, s in enumerate(silhouettes):
        if s >= best_silhouette - leeway:
            best_num = nums_to_try[idx]
    return best_num


# Check if there are overlaps and we should use cNMF rather than hierarchical in hybrid
def _check_overlap(X, clusters, cNMF_thresh, iters=20):
    X = X - 2
    cluster_no = len(set(clusters))

    # Only applicable if more than 3 clusters, so just return something that will return false
    if cluster_no < 3:
        return False

    # The representative vector is just a sum of the vectors from each cluster
    representative_vectors = np.zeros((cluster_no, len(X)))

    gene_counts = np.zeros(cluster_no)

    # Add each vector to the representative vector
    for cluster_i in range(1, cluster_no + 1):
        for cc, k in enumerate(clusters):
            if k == cluster_i:
                representative_vectors[k - 1] += X.iloc[:, cc].values
                gene_counts[k - 1] += 1

    rng = np.random.default_rng(42)

    cluster_means = representative_vectors / gene_counts[:, None]

    errors = []

    for i in range(representative_vectors.shape[0]):
        rel_test_error = 0
        for repeat in range(iters):
            target = cluster_means[i]
            others = np.delete(cluster_means, i, axis=0)

            # Only use coordinates where average gene count is > 3
            valid_idx = np.where((target > 3))[0]

            if len(valid_idx) < 2:
                continue  # not enough points to split

            train_idx = rng.choice(valid_idx, size=max(1, len(valid_idx) // 2), replace=False)

            test_mask = np.isin(valid_idx, train_idx, invert=True)
            test_idx = valid_idx[test_mask]

            if len(test_idx) == 0:
                continue

            # Fit on train coordinates only
            coeffs, *_ = np.linalg.lstsq(others[:, train_idx].T, target[train_idx], rcond=None)

            # Predict full vector
            prediction = coeffs @ others

            # Evaluate on held-out coordinates
            rel_test_error += np.linalg.norm(target[test_idx] - prediction[test_idx]) / np.sqrt(
                len(test_idx)
            )

            rel_test_error /= iters

        errors.append(rel_test_error)

    lowest_three = heapq.nsmallest(3, errors)

    return all(x < cNMF_thresh for x in lowest_three)
