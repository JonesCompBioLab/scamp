"""
Command-line tools for scAmp.
"""

from __future__ import annotations

import os
from pathlib import Path

from typing import Annotated, Literal, Union
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import scanpy as sc
import typer

from scamp import atac_cnv, vis, plotting, predict

scamp_app = typer.Typer(help="Tools for single-cell analysis of ecDNA.")

AnnDataFileArg = Annotated[
    str, typer.Option(help="File path to Anndata with copy-number data")
]
CopyNumberRangesDirArg = Annotated[
    str,
    typer.Option(
        help="File path to directory of RDSs of GRanges file of "
        "copy-number bins, if it's already computed."
    ),
]
CopyNumberFileArg = Annotated[
    str,
    typer.Argument(
        help="File path to anndata, tab/comma-delimited file, or MEX folder of copy number data"
    ),
]
FragDirArg = Annotated[
    str,
    typer.Option(
        help="Path to directory containing ATAC fragment files. If only one file, use --fragment-file"
    ),
]
ModelDirArg = Annotated[
    str, typer.Argument(help="Path to saved model directory.")
]
OutputDirArg = Annotated[str, typer.Argument(help="Directory of output files.")]
WhitelistFileArg = Annotated[
    str, typer.Option(help="File path to cellBC whitelist.")
]


@scamp_app.command(name="atac-cnv-legacy", help="Quantify single-cell copy-numbers.")
def quantify_copy_numbers_legacy(
    output_directory: OutputDirArg,
    copy_number_directory: CopyNumberRangesDirArg = None,
    fragment_directory: FragDirArg = None,
    whitelist_file: WhitelistFileArg = None,
    window_size: Annotated[
        int, typer.Option(help="Base pair width for genomic windows")
    ] = 3000000,
    step_size: Annotated[
        int,
        typer.Option(help="Step size from previous genomic window."),
    ] = 1000000,
    n_neighbors: Annotated[
        int,
        typer.Option(
            help="Number of genomic windows to compare against for normalization"
        ),
    ] = 200,
    reference_genome_name: Annotated[
        str,
        typer.Option(help="Reference genome name, to pair with a blacklist."),
    ] = "hg38",
):
    print('Warning: This is a legacy command and may be deprecated in future '
        'versions in favor of the python-based implementation `atac-cnv`.')

    binned_copy_number_script = (
        f"{os.path.dirname(__file__)}/scripts/scATAC_CNV.R"
    )
    gene_aggregation_script = (
        f"{os.path.dirname(__file__)}/scripts/aggregate_gene_copy_number.R"
    )

    # compute copy-numbers in genomic windows
    if fragment_directory:
        
        if whitelist_file is None:
            raise CLIError("If starting the copy-number pipeline from the " 
                            "beginning, please provide a whitelist file")
        print(
            f"Binning copy-numbers from {fragment_directory} in windows "
            f"of size {window_size}..."
        )
        os.system(
            f"Rscript {binned_copy_number_script} "
            f"{fragment_directory} {window_size} "
            f"{step_size} {n_neighbors} "
            f"{whitelist_file} {output_directory} "
            f"{reference_genome_name} {os.path.dirname(__file__)}"
        )

    # bin by gene
    if not copy_number_directory:
        copy_number_directory = output_directory

    print(f"Aggregating together copy-numbers across genes...")
    os.system(
        f"Rscript {gene_aggregation_script} "
        f"{copy_number_directory} {output_directory} {reference_genome_name} {os.path.dirname(__file__)}"
    )

@scamp_app.command(name="atac-cnv", help="Quantify single-cell copy-numbers.")
def quantify_copy_numbers(
    output_directory: Annotated[
        str,
        typer.Argument(help="Path to directory containing ATAC fragment files"),
    ] = None,
    fragment_directory: Annotated[
        str,
        typer.Option(
            help="Path to directory containing ATAC fragment files. If only one file, use --fragment-file"
        ),
    ] = None,
    fragment_file: Annotated[
        str, typer.Option(help="Path to ATAC fragment file")
    ] = None,
    sample_name: Annotated[
        str,
        typer.Option(
            help="If only using one file, use this instead of --fragment-file-key to specify sample name if file name format does not match sample name"
        ),
    ] = None,
    whitelist_file: WhitelistFileArg = None,
    pickle_dir: Annotated[
        str,
        typer.Option(
            help="Path to directory of pkl cell by window counts. pkl files will save here. If pkl files are "
            "present, skips the cell by window creation step"
        ),
    ] = None,
    window_size: Annotated[
        int, typer.Option(help="Base pair width for genomic windows")
    ] = 3000000,
    step_size: Annotated[
        int,
        typer.Option(help="Step size from previous genomic window."),
    ] = 1000000,
    n_neighbors: Annotated[
        int,
        typer.Option(
            help="Number of genomic windows to compare against for normalization"
        ),
    ] = 200,
    fragment_file_key: Annotated[
        str,
        typer.Argument(
            help="Path to tab separated two column list of matched file names and sample names"
        ),
    ] = None,
    cores_per_sample: Annotated[
        int, typer.Option(help="Number of cores to allocate per sample")
    ] = 16,
    max_workers: Annotated[
        int,
        typer.Option(
            help="Maximum number of workers (limit for memory bottlenecks. If not specified, becomes #cores/cores_per_sample)"
        ),
    ] = None,
    recreate_pkl: Annotated[
        bool,
        typer.Option(
            help="Rewrites pkl output file (use if edits were made to the sample with the same name and output location)"
        ),
    ] = True,
    reference_genome_name: Annotated[
        Literal["hg38", "hg19", "mm10", "mm39"],
        typer.Option(
            help="Reference genome name, to pair with a blacklist."
        ),
    ] = "hg38"
):
    if fragment_file is None and fragment_directory is None:
        print("Error: requires copy-numbers-file or copy-numbers-folder")

    import multiprocessing

    TOTAL_CORES = multiprocessing.cpu_count()

    if fragment_directory is not None:

        cores_per_sample = min(cores_per_sample, TOTAL_CORES)
        print(f"Using {cores_per_sample} cores for each sample")
        if max_workers is None:
            max_workers = max(1, TOTAL_CORES // cores_per_sample)
            print(f"Maximum workers active: {max_workers}")
            os.environ["OMP_NUM_THREADS"] = str(cores_per_sample)
            os.environ["MKL_NUM_THREADS"] = str(cores_per_sample)
            os.environ["OPENBLAS_NUM_THREADS"] = str(cores_per_sample)

    else:
        cores_per_sample = TOTAL_CORES

    atac_cnv.sc_cnv_pipeline(
        fragment_file,
        sample_name,
        fragment_directory,
        cores_per_sample,
        max_workers,
        whitelist_file,
        window_size,
        step_size,
        n_neighbors,
        output_directory,
        pickle_dir,
        recreate_pkl,
        fragment_file_key,
        reference_genome_name
    )

@scamp_app.command(
    name="visualize", help="Visualize ecDNA results with cellxgene"
)
def visualize(
    copy_numbers_file: Annotated[
        str,
        typer.Argument(
            help="Path to copy number data or copy number MEX folder"
        ),
    ],
    expression_file: Annotated[
        str,
        typer.Argument(
            help="Path to the expression data or expression MEX folder"
        ),
    ],
    scamp_tsv: Annotated[str, typer.Argument(help="Scamp Predict tsv")],
    umap_name: Annotated[
        str, typer.Option(help="Name of UMAP obsm in expression anndata")
    ] = "X_umap",
    temp_folder: Annotated[
        str, typer.Option(help="Folder for temporary anndata and scamp csv")
    ] = "./temp",
    cn_threshold: Annotated[
        float,
        typer.Option(
            help="Threshold for copy number for visualizing ecDNA genes. Set to -1 to not use"
        ),
    ] = 5,
    cn_percentile_threshold: Annotated[
        float,
        typer.Option(
            help="Threshold for copy number percentile for visualization. Leave default to not use"
        ),
    ] = 100,
):

    # Where the files will go
    os.makedirs(temp_folder, exist_ok=True)

    # Parse copy number data
    if os.path.isdir(copy_numbers_file):
        # MEX format
        cn_adata = sc.read_10x_mtx(copy_numbers_file)
    else:
        # If copy number, convert to anndata first
        copy_numbers_ext = copy_numbers_file.split(".")[-1]
        if copy_numbers_ext == "h5ad":
            cn_adata = vis.read_adata(copy_numbers_file)
        else:
            cn_adata = vis.setup_copynumber(copy_numbers_file)

    # Parse expression data
    if os.path.isdir(expression_file):
        # MEX format
        exp_adata = sc.read_10x_mtx(exp_adata)
    else:
        expression_file_ext = expression_file.split(".")[-1]
        if expression_file_ext == "h5ad":
            exp_adata = vis.read_adata(expression_file)
        else:
            exp_adata = vis.setup_expression(expression_file, cn_adata)

    # Get full anndata
    vis.setup_anndata(
        cn_adata,
        scamp_tsv,
        temp_folder,
        cn_threshold,
        cn_percentile_threshold,
        umap_name,
        exp_adata,
    )

    # Run cellxgene
    os.system(
        f"cellxgene launch {temp_folder}/annotated_anndata.h5ad --gene-sets-file {temp_folder}/ecDNA_gene_set.csv --open"
    )


@scamp_app.command(name="predict", help="Predict ecDNA status.")
def predict_ecdna(
    output_dir: OutputDirArg,
    copy_numbers: Annotated[
        str, typer.Argument(help="MEX, csv, tsv, or anndata of copy numbers (cell by gene matrix)")
    ],
    model_file: ModelDirArg = str(Path(__file__).resolve().parent / "pretrained_models" / "scamp_model_1.0"),


    whitelist_file: WhitelistFileArg = None,
    log_dir: Annotated[
        str, typer.Argument(help="cNMF log dir")
    ] = "./temp",

    verbose: Annotated[
        bool, typer.Option("--verbose", help="Enable verbose output")
    ] = False,
    decision_rule: Annotated[
        float, typer.Option(help="NN likelihood decision rule.")
    ] = 0.5,
    # min_copy_number: Annotated[
    #     float, typer.Option(help="Minimum copy-number to consider.")
    # ] = 2.0,
    # max_percentile: Annotated[
    #     float, typer.Option(help="Maximum percentile to cap copy-numbers.")
    # ] = 99.0,
    filter_copy_number: Annotated[
         float, typer.Option(help="Drop genes whose mean copy-number is below this threshold.")
    ] = 2.5,
    predictions: Annotated[
            str, typer.Option(help="Path to old predictions file if recalculating clusters/species without predicting again.")
    ] = None,
    predict_method: Annotated[
            Literal["NN", "KNN", "GMM"],
            typer.Option(help="Method to use for predicting ecDNA. NN (old version), KNN, or GMM")
    ] = "GMM",
    deconvolution_method: Annotated[
            Literal["hier", "cNMF", "auto"],
            typer.Option(help="Method to use for clustering genes into species. hier (hierarchical clustering based), cNMF, or auto (automatically chooses method)")
    ] = "auto",
    one_thresh: Annotated[
            float, typer.Option(help="Threshold for choosing one distribution in GMM/KMM + GMM.")
    ] = 2.5,
    kneedle_coeff: Annotated[
            float, typer.Option(help="Higher means kneedle tends to fewer GMM distributions.")
    ] = 3,
    hier_ddist: Annotated[
            float, typer.Option(help="Distance for null nodes in hierarchical clustering.")
    ] = 1.3,
    min_weight: Annotated[
            float, typer.Option(help="Distributions with lower weight are ignored in GMM.")
    ] = 0.005,
    var_scale: Annotated[
            float, typer.Option(help="Scaled GMM parameter. Greater value decreases false positives.")
    ] = 0.5,
    k_mult: Annotated[
            float, typer.Option(help="Multiplier for K in KNN for predictions. If not set, defaults to k = 100")
    ] = None,
    cNMF_thresh: Annotated[
            float, typer.Option(help="Threshold for choosing cNMF when using combo")
    ] = 0.55,
    error_w: Annotated[
            float, typer.Option(help="Importance of error in deciding cNMF number of ecDNA species")
    ] = 0.25,
    score_cutoff: Annotated[
            float, typer.Option(help="Cutoff for deciding a gene belongs to a species of ecDNA")
    ] = 3
) -> None:
    

    os.makedirs(output_dir, exist_ok=True)


    out_log = predict.run_sample(copy_numbers, output_dir, predictions, model_file, whitelist_file, decision_rule, 
               filter_copy_number, predict_method, one_thresh, kneedle_coeff, 
               min_weight, var_scale, k_mult, deconvolution_method, hier_ddist, cNMF_thresh,
               error_w, score_cutoff, log_dir, verbose) 

    # predictions.to_csv(f"{output_dir}/model_predictions.tsv", sep='\t')
    # if not no_plot:
    #     plotting.plot_scamp_predictions_plotly(
    #         predictions,
    #         f"{output_dir}/ecDNA_predictions.html",
    #         title=f"scAmp predictions for {copy_numbers_file.split('/')[-1]}"
    #     )

    with open(f"{output_dir}/log.log", "w") as f:
        f.writelines(line + "\n" for line in out_log)
    print(f"Output written out to {output_dir}.")
