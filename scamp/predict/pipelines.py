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
    one_thresh = None,
    kneedle_coeff = None,
    min_weight = None,
    var_scale = None,
    k_mult = None,
    verbose = False
) :
    model = models.SCAMP.load(saved_model_directory)

    if whitelist_file:
        whitelist = pd.read_csv(whitelist_file, header=None).iloc[:,0].values
        counts_df = counts_df.loc[np.intersect1d(counts_df.index, whitelist)]

    result_rows = []

    # Run one of three prediction methods
    # Keep track of vectors to avoid redundancy in computation
    used_vectors = {}
    for gene in counts_df.columns:
        if np.mean(counts_df[gene]) > filter_copy_number :
            rounded = tuple(np.round(counts_df[gene].to_numpy(), 3))

            # Already has near identical vector in database
            if rounded in used_vectors :
                new_row = used_vectors[rounded].copy()
                out_log.append(f"Skipping {gene}, Using: {used_vectors[rounded]}")
                new_row["Gene"] = gene
                result_rows.append(new_row)
            else :
                if predict_method == "NN" :
                    prediction, proba = predict_methods.NN_ecDNA_predict(model, counts_df[gene], out_log, np.array([gene]), decision_rule, verbose)
                    result_rows.append({"Gene" : gene, "Proba" : proba, "Prediction" : prediction})
                    used_vectors[rounded] = {"Gene" : gene, "Proba" : proba, "Prediction" : prediction}
                elif predict_method == "GMM" :
                    prediction, mean, var, weight = predict_methods.GMM_ecDNA_predict(model, counts_df[gene], out_log, one_thresh=one_thresh, genes = np.array([gene]), 
                                                                                    decision_rule=decision_rule, kneedle_coeff=kneedle_coeff, var_scale = var_scale, 
                                                                                    min_weight = min_weight, verbose = verbose)
                    result_rows.append({"Gene" : gene, "DistVar" : var, "DistMean" : mean, "Weight" : weight, "Prediction" : prediction})
                    used_vectors[rounded] = {"Gene" : gene, "DistVar" : var, "DistMean" : mean, "Weight" : weight, "Prediction" : prediction}

                elif predict_method == "KNN" :
                    prediction, num_ecDNA_pos_cells = predict_methods.KNN_ecDNA_predict(model, counts_df[gene], out_log, np.array([gene]), one_thresh = one_thresh,
                                                                                kneedle_coeff=kneedle_coeff, k_mult = k_mult,
                                                                                var_scale = var_scale, decision_rule = decision_rule, verbose = verbose)
                    result_rows.append({"Gene" : gene, "Num PosCells" : num_ecDNA_pos_cells, "Prediction" : prediction})
                    used_vectors[rounded] = {"Gene" : gene, "Num PosCells" : num_ecDNA_pos_cells, "Prediction" : prediction}
                else :
                    print("Predict method must be one of KNN, GMM, NN. We recommend using KNN")
                    exit(0)

    # In case no genes are above threshold
    if len(result_rows) == 0 :
        print("No genes passed copy number threshold")
        exit(0)

    prediction_df = pd.DataFrame(result_rows)

    return prediction_df



def cluster(
    filename, predict_df, counts_df, deconvolution_method, hier_ddist, cNMF_thresh, error_w, score_cutoff, log_dir, out_log, verbose
) :
    # Get only the ecDNA positive genes
    predict_df['Species'] = "None"
    ecDNA_genes = predict_df[predict_df["Prediction"] == True]["Gene"].tolist()
    ecDNA_genes = list(dict.fromkeys(ecDNA_genes))

    # If no genes are ecDNA positive, just return
    if len(ecDNA_genes) == 0 :
        return predict_df
    
    counts_df_ecDNA = counts_df[ecDNA_genes]

    # Get rid of duplicate copy number vectors
    df_unique, column_groups = species_deconvolution.remove_duplicates(counts_df_ecDNA)


    if deconvolution_method == "hier" :
        species_to_gene, cell_by_ecDNA = species_deconvolution.hier_deconvolution(df_unique, hier_ddist)
    elif deconvolution_method == "auto" :
        species_to_gene, cell_by_ecDNA = species_deconvolution.combo_deconvolution(df_unique, out_log, cNMF_thresh, filename, error_w = error_w, 
                                                                                   score_cutoff = score_cutoff, log_dir = log_dir, hier_ddist = hier_ddist, verbose = verbose)
    elif deconvolution_method == "cNMF" :
        species_to_gene, cell_by_ecDNA = species_deconvolution.cNMF_deconvolution(df_unique, filename, out_log, error_w = error_w, score_cutoff= score_cutoff,
                                                                                  log_dir = log_dir, hier_ddist = hier_ddist, verbose = verbose)
    else :
        print ("Cluster method must be one of hier, combo, or cNMF. We recommend using combo")
        exit(0)

    # Add labels back
    for species, genes in species_to_gene.items() :
        for gene in genes :
            for mapped_gene in column_groups[gene] :
                mask = predict_df["Gene"] == mapped_gene

                current = predict_df.loc[mask, "Species"].iloc[0]

                if current == "None":
                    predict_df.loc[mask, "Species"] = str(species)
                else:
                    predict_df.loc[mask, "Species"] = current + ',' + str(species)
    return predict_df, cell_by_ecDNA
    
    


# Can skip the first step by providing a predictions path (otherwise set it to none)
def run_sample(file, output_dir, predictions, model_file, whitelist_file, decision_rule, 
               filter_copy_number, predict_method, one_thresh, kneedle_coeff, 
               min_weight, var_scale, k_mult, deconvolution_method, hier_ddist, cNMF_thresh,
               error_w, score_cutoff, log_dir, verbose) :
    out_log = [f"Settings\nModel path: {model_file}\nDecision rule: {decision_rule}\nPrediction method: {predict_method}\nDeconvolution method: {deconvolution_method}"]
    out_log.append(f"Copy number filter: {filter_copy_number}")

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
        counts_df = io.read_copy_numbers_file(file)
    elif mode == "MEX" :
        counts_df = io.read_mex_file(file)
    else :
        counts_df = io.read_anndata_file(file)

    # Allow cutting in the middle of pipeline with old predictions file
    if predictions is None :
        prediction_df = predict(out_log, counts_df, model_file, whitelist_file, decision_rule, filter_copy_number, predict_method,
                                one_thresh, kneedle_coeff, min_weight, var_scale, k_mult, verbose)
    else :
        prediction_df = pd.read_csv(predictions, sep = '\t')

    # Output early in case something breaks
    filename = Path(file).stem
    prediction_df.to_csv(f"{output_dir}/ecDNA_preds_{filename}.tsv", sep='\t')


    prediction_df, cell_by_eCDNA = cluster(filename, prediction_df, counts_df, deconvolution_method, hier_ddist, cNMF_thresh, error_w, score_cutoff, log_dir, out_log, verbose)

    # Output predictions and visualizations
    prediction_df.to_csv(f"{output_dir}/ecDNA_preds_{filename}.tsv", sep='\t', index = False)
    cell_by_eCDNA.to_csv(f"{output_dir}/cell_by_ecDNA_{filename}.tsv", sep='\t', index = False)

    end = time.time()
    out_log.append(f"File {file} completed in {end - start:.2f} seconds")
    return out_log

