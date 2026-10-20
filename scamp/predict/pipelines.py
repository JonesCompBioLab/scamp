"""
A pipeline for predicting ecDNA status from single-cell copy-number
distributions.
"""

import warnings

warnings.filterwarnings("ignore", category=FutureWarning)

import os
import numpy as np
import pandas as pd
import scanpy as sc
import time
from pathlib import Path

from scamp import io
from scamp import models
from scamp.predict import predict_methods, species_deconvolution


def predict(
    out_log,
    counts_df,
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
        predict_method: One of "NN", "GMM" or "KNN".
        one_thresh: Threshold for choosing a single distribution in the GMM fit
            (GMM and KNN).
        kneedle_coeff: Higher values make kneedle pick fewer GMM distributions
            (GMM and KNN).
        min_weight: GMM distributions with less weight than this are ignored (GMM).
        var_scale: Distribution variances are scaled as var * weight ** var_scale.
            Larger values reduce false positives (GMM and KNN).
        k: Number of nearest neighbors for KNN.
        verbose: Whether to log per-gene details.

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
                if predict_method == "NN":
                    prediction, proba = predict_methods.NN_ecDNA_predict(
                        model, counts_df[gene], out_log, np.array([gene]), decision_rule, verbose
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

    return prediction_df


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
    # this is the usage hier_deconvolution and _hier_nnls give for a single profile.
    if df_unique.shape[1] == 1:
        species_to_gene = {"pred_ecDNA_1": list(df_unique.columns)}
        cell_by_ecDNA = (df_unique - 2).clip(lower=0)
        cell_by_ecDNA.columns = ["pred_ecDNA_1"]
    elif deconvolution_method == "hier":
        species_to_gene, cell_by_ecDNA = species_deconvolution.hier_deconvolution(
            df_unique, hier_ddist
        )
    elif deconvolution_method == "auto":
        species_to_gene, cell_by_ecDNA = species_deconvolution.combo_deconvolution(
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
    return predict_df, cell_by_ecDNA


# Can skip the first step by providing a predictions path (otherwise set it to none)
def run_sample(
    file,
    output_dir,
    predictions,
    model_file,
    whitelist_file,
    decision_rule,
    filter_copy_number,
    predict_method,
    one_thresh,
    kneedle_coeff,
    min_weight,
    var_scale,
    k,
    deconvolution_method,
    hier_ddist,
    cNMF_thresh,
    error_w,
    score_cutoff,
    log_dir,
    verbose,
):
    """Runs prediction and species deconvolution on one copy-number file.

    Writes `ecDNA_preds_<name>.tsv` and `cell_by_ecDNA_<name>.tsv` to `output_dir`,
    where `<name>` is the input file's stem. The predictions file is also written
    before deconvolution, so it is kept if deconvolution fails.

    Args:
        file: Path to a .h5ad, .csv or .tsv file, or a MEX directory (containing
            matrix.mtx* and barcodes.tsv*). Any other input is skipped and noted
            in the log.
        output_dir: Directory where the output TSVs are written.
        predictions: Path to an existing predictions TSV. If given, prediction is
            skipped and only deconvolution runs. Use None to predict.
        model_file: Path to a saved SCAMP model.
        whitelist_file: Headerless file with one cell barcode per line. Only the
            listed cells are kept. If None, all cells are used.
        decision_rule: SCAMP likelihood threshold for calling a gene ecDNA-positive.
        filter_copy_number: Drop genes whose mean copy-number is below this.
        predict_method: One of "NN", "GMM" or "KNN".
        one_thresh: Threshold for choosing a single distribution in the GMM fit
            (GMM and KNN).
        kneedle_coeff: Higher values make kneedle pick fewer GMM distributions
            (GMM and KNN).
        min_weight: GMM distributions with less weight than this are ignored (GMM).
        var_scale: Distribution variances are scaled as var * weight ** var_scale.
            Larger values reduce false positives (GMM and KNN).
        k: Number of nearest neighbors for KNN.
        deconvolution_method: One of "hier", "cNMF" or "auto" ("auto" chooses
            between the two).
        hier_ddist: Distance for null nodes in hierarchical clustering.
        cNMF_thresh: Threshold for choosing cNMF when the method is "auto".
        error_w: Weight of reconstruction error when choosing the number of cNMF species.
        score_cutoff: Score above which a gene is assigned to a species.
        log_dir: Directory for cNMF intermediate files.
        verbose: Whether to log extra details.

    Returns:
        The list of log lines for this run.
    """
    out_log = [
        f"Settings\nModel path: {model_file}\nDecision rule: {decision_rule}\nPrediction method: {predict_method}\nDeconvolution method: {deconvolution_method}"
    ]
    out_log.append(f"Copy number filter: {filter_copy_number}")

    # Detect extension
    out_log.append(f"Running {file}")
    p = Path(file)
    if os.path.isdir(p):
        has_matrix = any(p.glob("matrix.mtx*"))
        has_barcodes = any(p.glob("barcodes.tsv*"))
        if has_matrix and has_barcodes:
            mode = "MEX"
        else:
            out_log.append(f"{file} is non-MEX folder, skipping")
            return out_log

    else:
        copy_numbers_ext = file.split(".")[-1]
        if copy_numbers_ext == "h5ad":
            mode = "anndata"
        elif copy_numbers_ext == "csv" or copy_numbers_ext == "tsv":
            mode = "copynumber"
        else:
            out_log.append(f"{file} does not have extension h5ad, tsv, or csv. Skipping")
            return out_log

    out_log.append(f"Running {file}")
    start = time.time()
    # Call different wrapper for each prediction type
    if mode == "copynumber":
        counts_df = io.read_copy_numbers_file(file)
    elif mode == "MEX":
        counts_df = io.read_mex_file(file)
    else:
        counts_df = io.read_anndata_file(file)

    # Allow cutting in the middle of pipeline with old predictions file
    if predictions is None:
        prediction_df = predict(
            out_log,
            counts_df,
            model_file,
            whitelist_file,
            decision_rule,
            filter_copy_number,
            predict_method,
            one_thresh,
            kneedle_coeff,
            min_weight,
            var_scale,
            k,
            verbose,
        )
    else:
        prediction_df = pd.read_csv(predictions, sep="\t")

    # Output early in case something breaks
    filename = Path(file).stem
    prediction_df.to_csv(f"{output_dir}/ecDNA_preds_{filename}.tsv", sep="\t", index=False)

    prediction_df, cell_by_eCDNA = cluster(
        filename,
        prediction_df,
        counts_df,
        deconvolution_method,
        hier_ddist,
        cNMF_thresh,
        error_w,
        score_cutoff,
        log_dir,
        out_log,
        verbose,
    )

    # Output predictions and visualizations
    prediction_df.to_csv(f"{output_dir}/ecDNA_preds_{filename}.tsv", sep="\t", index=False)
    cell_by_eCDNA.to_csv(f"{output_dir}/cell_by_ecDNA_{filename}.tsv", sep="\t", index=False)

    end = time.time()
    out_log.append(f"File {file} completed in {end - start:.2f} seconds")
    return out_log
