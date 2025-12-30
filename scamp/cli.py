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
from scamp import cnv

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
    reference_fasta: Annotated[
        str, typer.Option(help="Path to reference fasta file (required for Python-native mode)")
    ] = "reference.fa",
    gene_annotation: Annotated[
        str, typer.Option(help="Path to gene annotation BED/GTF used for aggregation (required for Python-native mode)")
    ] = None,
    quiet: Annotated[
        bool, typer.Option(help="Suppress progress output")
    ] = False,
    use_pysam: Annotated[
        bool, typer.Option(help="Use pysam for faster chromosome loading (if available)")
    ] = True,
):

    # bin by gene
    if not copy_number_directory:
        copy_number_directory = output_directory

    if not reference_fasta:
        raise CLIError("`reference_fasta` is required for Python-native mode")
    if gene_annotation is None:
        raise CLIError("`gene_annotation` is required for Python-native mode")
    
    try:
        file_output = cnv.run_python_pipeline(
            fragment_directory,
            output_directory,
            reference_fasta,
            gene_annotation,
            window_size=window_size,
            step_size=step_size,
            n_neighbors=n_neighbors,
            quiet=quiet,
            use_pysam=use_pysam,
        )
        print(f"Aggregated CNAs wrote {len(file_output)} files.")
    except Exception as exc:
        raise CLIError(f"Error while running Python-native CNV pipeline: {exc}")


@scamp_app.command(name="predict", help="Predict ecDNA status.")
def predict_ecdna(
    output_dir: OutputDirArg,
    model_file: ModelDirArg,
    copy_numbers_file: CopyNumberFileArg = None,
    anndata_file: AnnDataFileArg = None,
    whitelist_file: WhitelistFileArg = None,
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
            filter_copy_number,
            whitelist_file
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