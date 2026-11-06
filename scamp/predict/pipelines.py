"""
A pipeline for predicting ecDNA status from single-cell copy-number
distributions.
"""

import warnings

warnings.filterwarnings("ignore", category=FutureWarning)

import numpy as np
import pandas as pd
import scanpy as sc

from scamp import io
from scamp import models
from scamp.predict import predict_methods, species_deconvolution


def predict(
    counts_df,
    out_log,
    saved_model_directory,
    whitelist_file,
    decision_rule,
    filter_copy_number,
    predict_method,
    one_thresh=1.15,
    kneedle_coeff=3,
    min_weight=0.005,
    var_scale=0.5,
    k=100,
    verbose=False,
    min_copy_number=2.0,
    max_percentile=99.9,
    trim_before_gmm=False,
    n_draws=1,
    single_dist_scamp=False,
    knn_min_mean=2.5,
    knn_min_var=10,
):
    """Predicts ecDNA status for each gene from its copy-number distribution across cells.

    Genes whose mean copy-number is at or below `filter_copy_number` are dropped.
    Genes with a near-identical copy-number vector (equal after rounding to 3
    decimals) reuse the first such gene's result rather than being predicted again.

    Args:
        out_log: List of log lines; progress messages are appended to it.
        counts_df: Cell-by-gene copy-number matrix (cells as index, genes as columns).
        saved_model_directory: Path to a saved SCAMP model.
        whitelist_file: Headerless file with one cell barcode per line. Only the
            listed cells are kept. If None, all cells are used.
        decision_rule: SCAMP likelihood threshold for calling a gene ecDNA-positive.
        filter_copy_number: Drop genes whose mean copy-number is below this.
        predict_method: One of "SCAMP", "GMM" or "KNN".
        one_thresh: Threshold for choosing a single distribution in the GMM fit
            (GMM and KNN).
        kneedle_coeff: Higher values make kneedle pick fewer GMM distributions
            (GMM and KNN).
        min_weight: GMM distributions with less weight than this are ignored (GMM).
        var_scale: Distribution variances are scaled as var * weight ** var_scale.
            Larger values reduce false positives (GMM and KNN).
        k: Number of nearest neighbors for KNN.
        verbose: Whether to log per-gene details.
        min_copy_number: Minimum of a copy-number distribution before prediction:
            copy-numbers at or below this are ignored when the model scores it.
        max_percentile: Copy-numbers above this percentile are dropped when the
            model scores a distribution.
        trim_before_gmm: Whether to also apply `min_copy_number` and
            `max_percentile` to the gene's cells before the GMM fit (GMM and KNN).
        n_draws: Draws of 1,000 sampled cells per distribution; the call uses the
            mean probability (GMM).
        single_dist_scamp: If only one distribution has weight above `min_weight`,
            score the gene's untrimmed cells as published scAmp does, with softmax
            >= `decision_rule` (GMM and KNN).
        knn_min_mean: A neighbourhood is scored only if its mean copy-number
            exceeds this (-inf = no gate) (KNN).
        knn_min_var: A neighbourhood is scored only if its copy-number variance
            exceeds this (-inf = no gate) (KNN).

    Returns:
        A DataFrame with one row per gene that passed the filter. It has `gene`,
        `mean` and `var` (the gene's copy-number mean and variance across cells),
        a boolean `prediction` column, and columns that depend on the method:
        `Probability` (NN); `Probability` and `Weight` (GMM); `NumberOfPositiveCells`
        (KNN).

    Raises:
        ValueError: If `predict_method` is not recognized, or no gene passes
            `filter_copy_number`.
    """
    model = models.SCAMP.load(saved_model_directory)

    if whitelist_file:
        whitelist = pd.read_csv(whitelist_file, header=None).iloc[:, 0].values
        counts_df = counts_df.loc[np.intersect1d(counts_df.index, whitelist)]

    result_rows = []

    # Run one of three prediction methods
    # Keep track of vectors to avoid redundancy in computation
    used_vectors = {}
    for gene in counts_df.columns:
        if np.mean(counts_df[gene]) > filter_copy_number:
            rounded = tuple(np.round(counts_df[gene].to_numpy(), 3))

            # Already has near identical vector in database
            if rounded in used_vectors:
                new_row = used_vectors[rounded].copy()
                out_log.append(f"Skipping {gene}, Using: {used_vectors[rounded]}")
                new_row["gene"] = gene
                result_rows.append(new_row)
            else:
                if predict_method == "SCAMP":
                    prediction, proba = predict_methods.SCAMP_ecDNA_predict(
                        model,
                        counts_df[gene],
                        out_log,
                        np.array([gene]),
                        decision_rule,
                        verbose,
                        min_copy_number=min_copy_number,
                        max_percentile=max_percentile,
                        filter_copy_number=filter_copy_number,
                    )
                    result_rows.append(
                        {"gene": gene, "Probability": proba, "prediction": prediction}
                    )
                    used_vectors[rounded] = {
                        "gene": gene,
                        "Probability": proba,
                        "prediction": prediction,
                    }
                elif predict_method == "GMM":
                    prediction, proba, _, _, weight = predict_methods.GMM_ecDNA_predict(
                        model,
                        counts_df[gene],
                        out_log,
                        one_thresh=one_thresh,
                        genes=np.array([gene]),
                        decision_rule=decision_rule,
                        kneedle_coeff=kneedle_coeff,
                        var_scale=var_scale,
                        min_weight=min_weight,
                        verbose=verbose,
                        min_copy_number=min_copy_number,
                        max_percentile=max_percentile,
                        filter_copy_number=filter_copy_number,
                        trim_before_gmm=trim_before_gmm,
                        n_draws=n_draws,
                        single_dist_scamp=single_dist_scamp,
                    )
                    result_rows.append(
                        {
                            "gene": gene,
                            "Probability": proba,
                            "Weight": weight,
                            "prediction": prediction,
                        }
                    )
                    used_vectors[rounded] = {
                        "gene": gene,
                        "Probability": proba,
                        "Weight": weight,
                        "prediction": prediction,
                    }

                elif predict_method == "KNN":
                    prediction, num_ecDNA_pos_cells = predict_methods.KNN_ecDNA_predict(
                        model,
                        counts_df[gene],
                        out_log,
                        np.array([gene]),
                        one_thresh=one_thresh,
                        kneedle_coeff=kneedle_coeff,
                        k=k,
                        var_scale=var_scale,
                        decision_rule=decision_rule,
                        verbose=verbose,
                        min_copy_number=min_copy_number,
                        max_percentile=max_percentile,
                        filter_copy_number=filter_copy_number,
                        trim_before_gmm=trim_before_gmm,
                        min_weight=min_weight,
                        knn_min_mean=knn_min_mean,
                        knn_min_var=knn_min_var,
                        single_dist_scamp=single_dist_scamp,
                    )
                    result_rows.append(
                        {
                            "gene": gene,
                            "NumberOfPositiveCells": num_ecDNA_pos_cells,
                            "prediction": prediction,
                        }
                    )
                    used_vectors[rounded] = {
                        "gene": gene,
                        "NumberOfPositiveCells": num_ecDNA_pos_cells,
                        "prediction": prediction,
                    }
                else:
                    raise ValueError(
                        f"Unknown predict method {predict_method!r}; "
                        "expected one of NN, KNN, GMM."
                    )

    # In case no genes are above threshold
    if len(result_rows) == 0:
        raise ValueError(
            f"No genes had a mean copy-number above {filter_copy_number}. "
            "Check the input file, or lower --filter-copy-number."
        )

    prediction_df = pd.DataFrame(result_rows)

    # Copy-number mean and variance of each gene across cells
    gene_counts = counts_df[prediction_df["gene"]]
    prediction_df.insert(1, "mean", gene_counts.mean().to_numpy())
    prediction_df.insert(2, "var", gene_counts.var(ddof=0).to_numpy())

    return prediction_df, out_log


def cluster(
    filename,
    predict_df,
    counts_df,
    deconvolution_method,
    hier_ddist,
    cNMF_thresh,
    error_w,
    score_cutoff,
    log_dir,
    out_log,
    verbose,
):
    """Groups the ecDNA-positive genes into ecDNA species.

    Genes with identical copy-number profiles are collapsed before deconvolution
    and expanded again afterward. If only one distinct profile remains, it is
    assigned to a single species, "pred_ecDNA_1".

    Args:
        filename: Sample name, used to label cNMF runs.
        predict_df: Predictions as returned by `predict`. Modified in place.
        counts_df: Cell-by-gene copy-number matrix used for the predictions.
        deconvolution_method: One of "hier", "cNMF" or "auto" ("auto" chooses
            between the two).
        hier_ddist: Distance for null nodes in hierarchical clustering.
        cNMF_thresh: Threshold for choosing cNMF when the method is "auto".
        error_w: Weight of reconstruction error when choosing the number of cNMF species.
        score_cutoff: Score above which a gene is assigned to a species.
        log_dir: Directory for cNMF intermediate files.
        out_log: List of log lines; progress messages are appended to it.
        verbose: Whether to log details of the deconvolution.

    Returns:
        A tuple `(predict_df, cell_by_ecDNA)`. `predict_df` gains a `Species`
        column: "None" for genes with no species, otherwise a comma-separated list
        of species names. `cell_by_ecDNA` is a cell-by-species DataFrame of
        estimated ecDNA amounts, which is empty (no columns) when no gene is ecDNA-positive.

    Raises:
        ValueError: If `deconvolution_method` is not recognized.
    """
    # Get only the ecDNA positive genes
    predict_df["Species"] = "None"
    ecDNA_genes = predict_df[predict_df["prediction"] == True]["gene"].tolist()
    ecDNA_genes = list(dict.fromkeys(ecDNA_genes))

    # If no genes are ecDNA positive, just return
    if len(ecDNA_genes) == 0:
        return predict_df, pd.DataFrame(index=counts_df.index)

    counts_df_ecDNA = counts_df[ecDNA_genes]

    # Get rid of duplicate copy number vectors
    df_unique, column_groups = species_deconvolution.remove_duplicates(counts_df_ecDNA)

    # One distinct copy-number profile is one species. The methods below need at least two
    # profiles to compare (np.corrcoef of a single column is a scalar), so assign it directly;
    # this is the usage hierarchical_deconvolution and _hier_nnls give for a single profile.
    if df_unique.shape[1] == 1:
        species_to_gene = {"pred_ecDNA_1": list(df_unique.columns)}
        cell_by_ecDNA = (df_unique - 2).clip(lower=0)
        cell_by_ecDNA.columns = ["pred_ecDNA_1"]
    elif deconvolution_method == "hier":
        species_to_gene, cell_by_ecDNA = species_deconvolution.hierarchical_deconvolution(
            df_unique, hier_ddist
        )
    elif deconvolution_method == "auto":
        species_to_gene, cell_by_ecDNA = species_deconvolution.hybrid_deconvolution(
            df_unique,
            out_log,
            cNMF_thresh,
            filename,
            error_w=error_w,
            score_cutoff=score_cutoff,
            log_dir=log_dir,
            hier_ddist=hier_ddist,
            verbose=verbose,
        )
    elif deconvolution_method == "cNMF":
        species_to_gene, cell_by_ecDNA = species_deconvolution.cNMF_deconvolution(
            df_unique,
            filename,
            out_log,
            error_w=error_w,
            score_cutoff=score_cutoff,
            log_dir=log_dir,
            hier_ddist=hier_ddist,
            verbose=verbose,
        )
    else:
        raise ValueError(
            f"Unknown deconvolution method {deconvolution_method!r}; "
            "expected one of hier, cNMF, auto."
        )

    # Add labels back
    for species, genes in species_to_gene.items():
        for gene in genes:
            for mapped_gene in column_groups[gene]:
                mask = predict_df["gene"] == mapped_gene

                current = predict_df.loc[mask, "Species"].iloc[0]

                if current == "None":
                    predict_df.loc[mask, "Species"] = str(species)
                else:
                    predict_df.loc[mask, "Species"] = current + "," + str(species)
    return predict_df, cell_by_ecDNA, out_log
