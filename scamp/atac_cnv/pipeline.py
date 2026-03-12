import os
from pathlib import Path
import pandas as pd
import numpy as np
import pyranges as pr
from tqdm import tqdm
import pickle

from .cnv_utils import *
import time
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed


'''
CNV from scATAC Pipeline

FRAGMENT_DIRECTORY: dir with fragment files (reads all that end with .tsv.gz)
cores_per_sample: number of cores to use on each sample (there is some parallel processing architecture per sample)
max_workers: number of samples to run in parallel
WHITELIST_FILE: whitelist
WINDOW_SIZE: window sizes for counts
STEP_SIZE: difference between windows
N_NEIGHBORS: number of nearest neghbors to average for GC bias
OUTPUT_DIRECTORY: output location
PKL_OUTPUT_DIRECTORY: intermediate pickle file location
recreate_pkl: whether or not to rewrite pickle file
FRAGMENT_FILE_KEY: a file that denotes the sample names for each fragment (if none, assumes correct file formats)
GENES_ANNO: gene annotation path
REFERENCE_BLACKLIST: blacklist file
'''
def sc_cnv_pipeline(FRAGMENT_FILE, FRAGMENT_DIRECTORY, cores_per_sample, max_workers, WHITELIST_FILE, WINDOW_SIZE, STEP_SIZE, N_NEIGHBORS, OUTPUT_DIRECTORY, 
                    PKL_OUTPUT_DIRECTORY, recreate_pkl, FRAGMENT_FILE_KEY = None, GENES_ANNO = '../reference/geneAnnohg38.tsv', 
                    REFERENCE_BLACKLIST = '../reference/hg38.blacklist.bed.gz') :
    os.makedirs(OUTPUT_DIRECTORY, exist_ok=True)

    # Download genome & read in genes
    script_dir = Path(__file__).resolve().parent

    gene_anno_path = script_dir / GENES_ANNO
    gene_anno_path = gene_anno_path.resolve()
    genome_path = script_dir / '../reference/.fa_genomes'
    genome = get_hg38_fasta(genome_path)
    genes_pr = get_gene_pr(gene_anno_path)


    # {fragment file : sample name}
    if FRAGMENT_DIRECTORY is not None :
        frag_dict = read_frag_files(FRAGMENT_DIRECTORY, FRAGMENT_FILE_KEY)

    # Get whitelist
    whitelists = get_whitelists(WHITELIST_FILE)


    # Create windows with blacklist
    print("Getting windows")
    blacklist_path = script_dir / REFERENCE_BLACKLIST
    blacklist_path = blacklist_path.resolve()
    windows = get_windows(genome, WINDOW_SIZE, STEP_SIZE, blacklist_path)
    print("Windows analyzed")

    total_time_start = time.time()

    if FRAGMENT_DIRECTORY is not None :
        with ProcessPoolExecutor(max_workers = max_workers) as executor :
            logs = []
            for frag_file, sample_name in frag_dict.items() :
                logs.append(
                    executor.submit(
                        run_sample, frag_file, sample_name,
                        OUTPUT_DIRECTORY, PKL_OUTPUT_DIRECTORY,
                        windows, whitelists, N_NEIGHBORS, genes_pr,
                        recreate_pkl, cores_per_sample
                    )
                )

            for log in as_completed(logs):
                result = log.result()
                if result is not None:
                    for line in result :
                        print(line)
    else :
        path = Path(FRAGMENT_FILE)
        if '-atac' in path.name :
            sample_name = path.name.replace("-atac_fragments.tsv.gz", "")
        elif '_atac' in path.name :
            sample_name = path.name.replace("_atac_fragments.tsv.gz", "")
        else :
            sample_name = path.name.replace("-fragments.tsv.gz", "")
        log = run_sample(FRAGMENT_FILE, sample_name,
                        OUTPUT_DIRECTORY, PKL_OUTPUT_DIRECTORY,
                        windows, whitelists, N_NEIGHBORS, genes_pr,
                        recreate_pkl, cores_per_sample)
        for i in log:
            print(i)
        
    total_time_end = time.time()
    print(f"Total time for all samples: {total_time_end - total_time_start:.2f} seconds")

    print("Done")

def run_sample(frag_file, sample_name, OUTPUT_DIRECTORY, PKL_OUTPUT_DIRECTORY, windows, whitelists, N_NEIGHBORS, genes_pr, recreate_pkl, cores_per_sample) :
    out_log = []
    out_log.append("")
    out_log.append(f"Processing {sample_name}")
    # Pickle file output
    if PKL_OUTPUT_DIRECTORY is None :
        PKL_OUTPUT_DIRECTORY = f"{OUTPUT_DIRECTORY}/pkl_files"
        os.makedirs(PKL_OUTPUT_DIRECTORY, exist_ok=True)

    if not os.path.exists(PKL_OUTPUT_DIRECTORY):
        os.makedirs(PKL_OUTPUT_DIRECTORY, exist_ok=True)

    # Search for pickle file if already exists
    pickle_out = f"{PKL_OUTPUT_DIRECTORY}/{sample_name}_windowstats.pkl"
    if Path(pickle_out).exists() and not recreate_pkl :
        out_log.append(f"Pickle file found: {pickle_out}")

        with open(pickle_out, "rb") as f:
            data_package = pickle.load(f)
    # Otherwise run pipeline
    else :
        out_log.append("Creating cell by window count matrix")
        start = time.time()
        data_package = run_aggregation(out_log, frag_file, sample_name, pickle_out, windows, whitelists, N_NEIGHBORS, bgdCN=2, MAKE_TEMP_SAVE=True)
        end = time.time()
        out_log.append(f"Cell by window runtime: {end - start:.2f} seconds")


    if data_package == None:
        out_log.append(f"Count failed for {sample_name}")
        return out_log
    
    start = time.time()


    # Widen by 1e5
    data_package['wmeta']["Start"] = data_package['wmeta']["Start"] - 100000
    data_package['wmeta']["End"] = data_package['wmeta']["End"] + 100000
    data_package['wmeta']["Start"] = data_package['wmeta']["Start"].clip(lower=0)

    out_log.append("Aggregating windows for genes")
    gene_matrix, gene_to_idx = aggregate_genes(genes_pr, data_package)
    col_names = [name for name, idx in sorted(gene_to_idx.items(), key=lambda x: x[1])]
    row_names = [name for name, idx in sorted(data_package['barcodeidx'].items(), key=lambda x: x[1])]

    cnv_df = pd.DataFrame(gene_matrix.T, index=row_names, columns=col_names)

    # Export as TSV
    out_log.append("Exporting")
    cnv_df.to_csv(f"{OUTPUT_DIRECTORY}/{sample_name}_cnv.tsv", sep="\t")

    end = time.time()
    out_log.append(f"Gene aggregation runtime {end - start:.2f} seconds")
    out_log.append("")
    return(out_log)