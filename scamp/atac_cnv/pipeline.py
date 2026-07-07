"""
Pipeline entry points for ATAC-derived copy-number calling.
"""
import os

from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from pathlib import Path
import time

import numpy as np
import pandas as pd
import pyranges as pr
import pickle
from tqdm import tqdm

from . import cnv_utils

def sc_cnv_pipeline(
    FRAGMENT_FILE,
    single_sample_name,
    FRAGMENT_DIRECTORY,
    cores_per_sample,
    max_workers,
    WHITELIST_FILE,
    WINDOW_SIZE,
    STEP_SIZE,
    N_NEIGHBORS,
    OUTPUT_DIRECTORY,
    PKL_OUTPUT_DIRECTORY,
    recreate_pkl,
    FRAGMENT_FILE_KEY=None,
    ASSEMBLY="hg38"
):
    """Run the scATAC copy-number pipeline.

    Reads fragment files, creates or loads genome windows, computes window-level
    copy-number estimates, aggregates copy-number values to genes, and writes
    one CNV TSV per sample.

    Args:
        FRAGMENT_FILE: Path to a single fragments file when not using a
            fragment directory.
        single_sample_name: Optional sample name for a single fragments file.
        FRAGMENT_DIRECTORY: Directory with fragment files ending in .tsv.gz.
        cores_per_sample: Number of cores to use within each sample.
        max_workers: Number of samples to process in parallel.
        WHITELIST_FILE: Optional barcode whitelist file.
        WINDOW_SIZE: Width of each counting window in bases.
        STEP_SIZE: Distance in bases between consecutive windows.
        N_NEIGHBORS: Number of GC-nearest windows to average for GC bias.
        OUTPUT_DIRECTORY: Directory where CNV TSV files are written.
        PKL_OUTPUT_DIRECTORY: Directory where intermediate pickle files are
            written. If None, a pkl_files directory is created under the output
            directory.
        recreate_pkl: Whether to recreate existing intermediate pickle files.
        FRAGMENT_FILE_KEY: Optional file mapping fragment files to sample names.
        GENES_ANNO: Path to gene annotation file, relative to this module.
        REFERENCE_BLACKLIST: Path to blacklist BED file, relative to this
            module.
    """
    os.makedirs(OUTPUT_DIRECTORY, exist_ok=True)

    # Download genome & read in genes
    script_dir = Path(__file__).resolve().parent

    gene_anno_path = script_dir / f"../reference/geneAnno{ASSEMBLY}.tsv"
    gene_anno_path = gene_anno_path.resolve()
    genome_path = script_dir / f"../reference/.fa_genomes_{ASSEMBLY}"
    genome = cnv_utils.get_assembly_fasta(ASSEMBLY, genome_path)
    genes_pr = cnv_utils.get_gene_pr(gene_anno_path)

    # {fragment file : sample name}
    if FRAGMENT_DIRECTORY is not None:
        frag_dict = cnv_utils.read_frag_files(FRAGMENT_DIRECTORY, FRAGMENT_FILE_KEY)

    # Get whitelist
    whitelists = cnv_utils.get_whitelists(WHITELIST_FILE)

    # Create windows with blacklist
    print("Getting windows")
    blacklist_path = script_dir / f"../reference/{ASSEMBLY}-blacklist.bed.gz"
    blacklist_path = blacklist_path.resolve()
    windows = cnv_utils.get_windows(genome, WINDOW_SIZE, STEP_SIZE, blacklist_path)
    print("Windows analyzed")

    total_time_start = time.time()

    if FRAGMENT_DIRECTORY is not None:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            logs = []
            for frag_file, sample_name in frag_dict.items():
                logs.append(
                    executor.submit(
                        run_sample,
                        frag_file,
                        sample_name,
                        OUTPUT_DIRECTORY,
                        PKL_OUTPUT_DIRECTORY,
                        windows,
                        whitelists,
                        N_NEIGHBORS,
                        genes_pr,
                        recreate_pkl,
                        cores_per_sample,
                    )
                )

            for log in as_completed(logs):
                result = log.result()
                if result is not None:
                    for line in result:
                        print(line)
    else:
        # Get sample names
        path = Path(FRAGMENT_FILE)
        if single_sample_name is None:
            if "-atac" in path.name:
                sample_name = path.name.replace("-atac_fragments.tsv.gz", "")
            elif "_atac" in path.name:
                sample_name = path.name.replace("_atac_fragments.tsv.gz", "")
            else:
                sample_name = path.name.replace("-fragments.tsv.gz", "")
        else:
            sample_name = single_sample_name
        log = run_sample(
            FRAGMENT_FILE,
            sample_name,
            OUTPUT_DIRECTORY,
            PKL_OUTPUT_DIRECTORY,
            windows,
            WINDOW_SIZE,
            whitelists,
            N_NEIGHBORS,
            genes_pr,
            recreate_pkl,
            cores_per_sample,
        )
        for i in log:
            print(i)

    total_time_end = time.time()
    print(
        f"Total time for all samples: {total_time_end - total_time_start:.2f} seconds"
    )

    print("Done")


def run_sample(
    frag_file,
    sample_name,
    OUTPUT_DIRECTORY,
    PKL_OUTPUT_DIRECTORY,
    windows,
    window_size,
    whitelists,
    N_NEIGHBORS,
    genes_pr,
    recreate_pkl,
    cores_per_sample,
):
    """Run CNV aggregation and gene export for one sample.

    Args:
        frag_file: Path to a gzipped fragments TSV file.
        sample_name: Sample name used in output filenames and whitelist lookup.
        OUTPUT_DIRECTORY: Directory where the CNV TSV file is written.
        PKL_OUTPUT_DIRECTORY: Directory where intermediate pickle files are
            written. If None, a pkl_files directory is created under the output
            directory.
        windows: Dataframe of annotated windows.
        whitelists: Optional dictionary mapping sample names to barcodes.
        N_NEIGHBORS: Number of GC-nearest windows to average for GC bias.
        genes_pr: PyRanges object containing gene annotations.
        recreate_pkl: Whether to recreate existing intermediate pickle files.
        cores_per_sample: Number of cores to use within the sample.
    """
    out_log = []
    out_log.append("")
    out_log.append(f"Processing {sample_name}")
    # Pickle file output
    if PKL_OUTPUT_DIRECTORY is None:
        PKL_OUTPUT_DIRECTORY = f"{OUTPUT_DIRECTORY}/pkl_files"
        os.makedirs(PKL_OUTPUT_DIRECTORY, exist_ok=True)

    if not os.path.exists(PKL_OUTPUT_DIRECTORY):
        os.makedirs(PKL_OUTPUT_DIRECTORY, exist_ok=True)

    # Search for pickle file if already exists
    pickle_out = f"{PKL_OUTPUT_DIRECTORY}/{sample_name}_windowstats.pkl"
    if Path(pickle_out).exists() and not recreate_pkl:
        out_log.append(f"Pickle file found: {pickle_out}")

        with open(pickle_out, "rb") as f:
            data_package = pickle.load(f)
            
    # Otherwise run pipeline
    else:
        out_log.append("Creating cell by window count matrix")
        start = time.time()
        data_package = cnv_utils.run_aggregation(
            out_log,
            frag_file,
            sample_name,
            pickle_out,
            windows,
            window_size,
            whitelists,
            N_NEIGHBORS,
            bgdCN=2,
            MAKE_TEMP_SAVE=True,
        )
        end = time.time()
        out_log.append(f"Cell by window runtime: {end - start:.2f} seconds")

    if data_package == None:
        out_log.append(f"Count failed for {sample_name}")
        return out_log

    start = time.time()

    # Widen by 1e5
    data_package["wmeta"]["Start"] = data_package["wmeta"]["Start"] - 100000
    data_package["wmeta"]["End"] = data_package["wmeta"]["End"] + 100000
    data_package["wmeta"]["Start"] = data_package["wmeta"]["Start"].clip(
        lower=0
    )

    out_log.append("Aggregating windows for genes")
    gene_matrix, gene_to_idx = cnv_utils.aggregate_genes(genes_pr, data_package)
    col_names = [
        name for name, idx in sorted(gene_to_idx.items(), key=lambda x: x[1])
    ]
    row_names = [
        name
        for name, idx in sorted(
            data_package["barcodeidx"].items(), key=lambda x: x[1]
        )
    ]

    cnv_df = pd.DataFrame(gene_matrix.T, index=row_names, columns=col_names)

    # Export as TSV
    out_log.append("Exporting")
    cnv_df.to_csv(f"{OUTPUT_DIRECTORY}/{sample_name}_cnv.tsv", sep="\t")

    end = time.time()
    out_log.append(f"Gene aggregation runtime {end - start:.2f} seconds")
    out_log.append("")
    return out_log
