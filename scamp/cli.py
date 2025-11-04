"""
Command-line tools for scAmp.
"""

from __future__ import annotations

import os
import pathlib
import pickle
from typing import Annotated, Union

import typer

from scamp import io
from scamp.mixins import CLIError
from scamp import predict
from scamp import plotting
from scamp import vis

scamp_app = typer.Typer(help="Tools for single-cell analysis of ecDNA.")

AnnDataFileArg = Annotated[
    str, typer.Option(help="File path to Anndata with copy-number data")
]
CopyNumberRangesDirArg = Annotated[
    str, typer.Option(help="File path to directory of RDSs of GRanges file of "
                      "copy-number bins, if it's already computed.")
]
CopyNumberFileArg = Annotated[
    str, typer.Option(help="File path to tab-delimited file of copy numbers")
]
FragDirArg = Annotated[
    str,
    typer.Option(help="Path to directory containing ATAC fragment files."),
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
                            "beginning, please provide a whitelist fiel")
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
        f"{copy_number_directory} {output_directory} {reference_genome_name}"
    )

@scamp_app.command(name="visualize", help="Visualize ecDNA results with cellxgene")
def visualize(
    copy_numbers_file: CopyNumberFileArg = None,
    anndata_file: AnnDataFileArg = None,
    mode: Annotated[
        str, typer.Option(help="Mode: anndata copynumber")
    ] = "copynumber",
    scamp_tsv: Annotated[
        str, typer.Option(help="Scamp Predict tsv")
    ] = ...,
    umap_name: Annotated[
        str, typer.Option(help='Name of UMAP obsm in anndata')
    ] = "X_umap",
    expression_file: Annotated[
        str, typer.Option(help="If you want to visualize with ")
    ] = None,
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
    # Where the files will go
    os.makedirs(temp_folder, exist_ok=True)

    # If copy number, convert to anndata first
    if mode == "anndata" :
        adata = vis.read_adata(anndata_file)
        vis.setup_anndata(adata, scamp_tsv, temp_folder, cn_threshold, cn_percentile_threshold, umap_name, expression_file)
    else :
        vis.setup_copynumber(copy_numbers_file, scamp_tsv, temp_folder, cn_threshold, cn_percentile_threshold, umap_name, expression_file) 

    # Run cellxgene   
    os.system(f"cellxgene launch {temp_folder}/annotated_anndata.h5ad --gene-sets-file {temp_folder}/ecDNA_gene_set.csv --open")


    

@scamp_app.command(name="predict", help="Predict ecDNA status.")
def predict_ecdna(
    output_dir: OutputDirArg,
    model_file: ModelDirArg,
    copy_numbers_file: CopyNumberFileArg = None,
    anndata_file: AnnDataFileArg = None,
    mode: Annotated[
        str, typer.Option(help="Mode: (currently only offering `copynumber`)")
    ] = "copynumber",
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
    ] = False
) -> None:

    if (copy_numbers_file is None) and (anndata_file is None):
        raise CLIError("Specify one of copy numbers file anndata file.")

    if mode == "copynumber":
        predictions = predict.predict_ecdna_from_copy_number(
            copy_numbers_file,
            model_file,
            decision_rule,
            min_copy_number,
            max_percentile,
            filter_copy_number
        )

        os.makedirs(output_dir)

        predictions.to_csv(f"{output_dir}/model_predictions.tsv", sep='\t')
        if not no_plot:
            plotting.plot_scamp_predictions_plotly(
                predictions,
                f"{output_dir}/ecDNA_predictions.html",
                title=f"scAmp predictions for {copy_numbers_file.split('/')[-1]}"
            )

        print(f"Output written out to {output_dir}.")

