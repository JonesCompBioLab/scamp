"""
Command-line tools for scAmp.
"""
from __future__ import annotations

import os

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
from typing import Annotated, Union
import typer

scamp_app = typer.Typer(help="Tools for single-cell analysis of ecDNA.")

AnnDataFileArg = Annotated[
    str, typer.Option(help="File path to Anndata with copy-number data")
]
CopyNumberRangesDirArg = Annotated[
    str, typer.Option(help="File path to directory of RDSs of GRanges file of "
                      "copy-number bins, if it's already computed.")
]
CopyNumberFileArg = Annotated[
    str, typer.Argument(help="File path to anndata, tab/comma-delimited file, or MEX folder of copy number data")
]
FragDirArg = Annotated[
    str,
    typer.Option(help="Path to directory containing ATAC fragment files. If only one file, use --fragment-file"),
]
ModelDirArg = Annotated[
    str, typer.Argument(help="Path to saved model directory.")
]
OutputDirArg = Annotated[str, typer.Argument(help="Directory of output files.")]
WhitelistFileArg = Annotated[
    str, typer.Option(help="File path to cellBC whitelist.")
]


@scamp_app.command(name="atac-cnv", help="Quantify single-cell copy-numbers.")
def quantify_copy_numbers(
    output_directory: Annotated[
            str, typer.Argument(help="Path to directory containing ATAC fragment files")
    ] = None,     
    fragment_directory: Annotated[
        str,
        typer.Option(help="Path to directory containing ATAC fragment files. If only one file, use --fragment-file"),
    ] = None,
    fragment_file: Annotated[
            str, typer.Option(help="Path to ATAC fragment file")
    ] = None,
    sample_name: Annotated[
            str, typer.Option(help="If only using one file, use this instead of --fragment-file-key to specify sample name if file name format does not match sample name")
    ] = None,
    whitelist_file: WhitelistFileArg = None,

    pickle_dir: Annotated[
        str, typer.Option(help="Path to directory of pkl cell by window counts. pkl files will save here. If pkl files are " \
        "present, skips the cell by window creation step")
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
        str, typer.Argument(help="Path to tab separated two column list of matched file names and sample names")
    ] = None,
    cores_per_sample: Annotated[
        int, typer.Option(help="Number of cores to allocate per sample")
    ] = 16,
    max_workers: Annotated[
        int, typer.Option(help="Maximum number of workers (limit for memory bottlenecks. If not specified, becomes #cores/cores_per_sample)")
    ] = None,
    recreate_pkl: Annotated[
        bool,
        typer.Option(help="Rewrites pkl output file (use if edits were made to the sample with the same name and output location)"),
    ] = False,
    reference_genome_name: Annotated[
        str,
        typer.Option(help="Reference genome name, to pair with a blacklist. Currently only supports hg38"),
    ] = "hg38"
):
    if fragment_file is None and fragment_directory is None :
        print("Error: requires copy-numbers-file or copy-numbers-folder")

    import multiprocessing

    TOTAL_CORES = multiprocessing.cpu_count()

    if fragment_directory is not None :

        cores_per_sample = min(cores_per_sample, TOTAL_CORES)
        print(f"Using {cores_per_sample} cores for each sample")
        if max_workers is None :
            max_workers = max(1, TOTAL_CORES // cores_per_sample)
            print(f"Maximum workers active: {max_workers}")
            os.environ["OMP_NUM_THREADS"] = str(cores_per_sample)
            os.environ["MKL_NUM_THREADS"] = str(cores_per_sample)
            os.environ["OPENBLAS_NUM_THREADS"] = str(cores_per_sample)
        
    else :
        cores_per_sample = TOTAL_CORES

    from scamp import atac_cnv

    atac_cnv.sc_cnv_pipeline(fragment_file, sample_name, fragment_directory, cores_per_sample, max_workers, whitelist_file, window_size, step_size, 
                             n_neighbors, output_directory, pickle_dir, recreate_pkl, fragment_file_key, GENES_ANNO = '../reference/geneAnnohg38.tsv', 
                             REFERENCE_BLACKLIST = '../reference/hg38.blacklist.bed.gz')



@scamp_app.command(name="visualize", help="Visualize ecDNA results with cellxgene")
def visualize(
    copy_numbers_file: Annotated[
        str, typer.Argument(help="Path to copy number data or copy number MEX folder")
    ],
    expression_file: Annotated[
        str, typer.Argument(help="Path to the expression data or expression MEX folder")
    ],
    scamp_tsv: Annotated[
        str, typer.Argument(help="Scamp Predict tsv")
    ],
    umap_name: Annotated[
        str, typer.Option(help='Name of UMAP obsm in expression anndata')
    ] = "X_umap",
    temp_folder: Annotated[
        str, typer.Option(help="Folder for temporary anndata and scamp csv")
    ] = "./temp",
    cn_threshold: Annotated[
        float, typer.Option(help='Threshold for copy number for visualizing ecDNA genes. Set to -1 to not use')
    ] = 5,
    cn_percentile_threshold: Annotated[
        float, typer.Option(help='Threshold for copy number percentile for visualization. Leave default to not use')
    ] = 100


) :
    from scamp import vis
    import scanpy as sc

    # Where the files will go
    os.makedirs(temp_folder, exist_ok=True)

    # Parse copy number data
    if os.path.isdir(copy_numbers_file) :
        # MEX format
        cn_adata = sc.read_10x_mtx(copy_numbers_file)
    else :
        # If copy number, convert to anndata first
        copy_numbers_ext = copy_numbers_file.split('.')[-1]
        if copy_numbers_ext == "h5ad" :
            cn_adata = vis.read_adata(copy_numbers_file)
        else :
            cn_adata = vis.setup_copynumber(copy_numbers_file) 
    
    # Parse expression data
    if os.path.isdir(expression_file) :
        # MEX format
        exp_adata = sc.read_10x_mtx(exp_adata)
    else :
        expression_file_ext = expression_file.split('.')[-1]
        if expression_file_ext == "h5ad" :
            exp_adata = vis.read_adata(expression_file)
        else :
            exp_adata = vis.setup_expression(expression_file, cn_adata)

    # Get full anndata
    vis.setup_anndata(cn_adata, scamp_tsv, temp_folder, cn_threshold, cn_percentile_threshold, umap_name, exp_adata)

    # Run cellxgene   
    os.system(f"cellxgene launch {temp_folder}/annotated_anndata.h5ad --gene-sets-file {temp_folder}/ecDNA_gene_set.csv --open")


    

@scamp_app.command(name="predict", help="Predict ecDNA status.")
def predict_ecdna(
    output_dir: OutputDirArg,
    model_file: ModelDirArg,
    copy_numbers_file: Annotated[
        str, typer.Option(help="File path to anndata, tab/comma-delimited file, or MEX folder of copy number data")
    ] = None,
    copy_numbers_folder: Annotated[
        str, typer.Option(help="Folder path containing anndatas, tab/comma-delimited files, or MEX folders of copy number data")
    ] = None,
    whitelist_file: WhitelistFileArg = None,
    decision_rule: Annotated[
        float, typer.Option(help="Likelihood decision rule.")
    ] = 0.5,
    min_copy_number: Annotated[
        float, typer.Option(help="Minimum copy-number to consider.")
    ] = 2.0,
    max_percentile: Annotated[
        float, typer.Option(help="Maximum percentile to cap copy-numbers.")
    ] = 99.0,
    filter_copy_number: Annotated[
         float, typer.Option(help="Drop genes whose mean copy-number is below this threshold.")
    ] = 2.5,
    no_plot: Annotated[
        float, typer.Option(help="Suppress plotting functionality.")
    ] = False,
    cluster_distance_threshold: Annotated[
        float, typer.Option(help="Distance threshold for hierarchical clustering.")
    ] = 0.4,
    cores_per_sample: Annotated[
        int, typer.Option(help="Number of cores to allocate per sample")
    ] = 16,
    max_workers: Annotated[
        int, typer.Option(help="Maximum number of workers (limit for memory bottlenecks. If not specified, becomes #cores/cores_per_sample)")
    ] = None
) -> None:
    

    if copy_numbers_file is None and copy_numbers_folder is None :
        print("Error: requires copy-numbers-file or copy-numbers-folder")

    import multiprocessing

    TOTAL_CORES = multiprocessing.cpu_count()

    if copy_numbers_folder is not None :
        cores_per_sample = min(cores_per_sample, TOTAL_CORES)
        print(f"Using {cores_per_sample} cores for each sample")
        if max_workers is None :
            max_workers = max(1, TOTAL_CORES // cores_per_sample)
        print(f"Maximum workers active: {max_workers}")
        os.environ["OMP_NUM_THREADS"] = str(cores_per_sample)
        os.environ["MKL_NUM_THREADS"] = str(cores_per_sample)
        os.environ["OPENBLAS_NUM_THREADS"] = str(cores_per_sample)
    else :
        max_workers = 1
        cores_per_sample = TOTAL_CORES

    from scamp import predict
    from pathlib import Path
    from concurrent.futures import ProcessPoolExecutor, as_completed
    import time
    
    if copy_numbers_folder is None :
        # Just one file if only one provided
        files_list = []
        files_list.append(copy_numbers_file)
    else :
        data_dir = Path(copy_numbers_folder).resolve()
        files_list = [str(f) for f in data_dir.glob("*")]

    # Parallelization
    total_time_start = time.time()
    with ProcessPoolExecutor(max_workers = max_workers) as executor :
        logs = []
        for file in files_list :
            logs.append(
                executor.submit(
                    predict.run_sample, file, output_dir, model_file, whitelist_file, decision_rule, min_copy_number, max_percentile, 
                                        filter_copy_number, cluster_distance_threshold, no_plot)
                )
            
        for log in as_completed(logs):
            result = log.result()
            if result is not None:
                for line in result :
                    print(line)
    
    total_time_end = time.time()
    print(f"Total time for all samples: {total_time_end - total_time_start:.2f} seconds")

    print("Done")