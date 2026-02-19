from pathlib import Path
import pandas as pd
import numpy as np
import pyranges as pr
import os
from tqdm import tqdm
import pickle

from .cnv_utils import *
import time

# # TODO: currently, this whitelist is for one sample, typically the whitelist has samplename#cell, so try to resolve that!
# WHITELIST_FILE = "ML499M1-S1_whitelist.txt"
# WINDOW_SIZE = 3000000
# STEP_SIZE = 1000000
# N_NEIGHBORS = 200

'''
CNV from scATAC Pipeline

FRAGMENT_DIRECTORY: dir with fragment files (reads all that end with .tsv.gz)
WHITELIST_FILE: whitelist
WINDOW_SIZE: window sizes for counts
STEP_SIZE: difference between windows
N_NEIGHBORS: number of nearest neghbors to average for GC bias
OUTPUT_DIRECTORY: output location
PKL_OUTPUT_DIRECTORY: intermediate pickle file location
FRAGMENT_FILE_KEY: a file that denotes the sample names for each fragment (if none, assumes correct file formats)
GENES_ANNO: gene annotation path
REFERENCE_BLACKLIST: blacklist file
'''
def sc_cnv_pipeline(FRAGMENT_DIRECTORY, WHITELIST_FILE, WINDOW_SIZE, STEP_SIZE, N_NEIGHBORS, OUTPUT_DIRECTORY, 
                    PKL_OUTPUT_DIRECTORY, FRAGMENT_FILE_KEY = None, GENES_ANNO = '../reference/geneAnnohg38.tsv', 
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
    frag_dict = read_frag_files(FRAGMENT_DIRECTORY, FRAGMENT_FILE_KEY)

    # Get whitelist
    whitelists = get_whitelists(WHITELIST_FILE)


    # Create windows with blacklist
    print("Creating windows")
    blacklist_path = script_dir / REFERENCE_BLACKLIST
    blacklist_path = blacklist_path.resolve()
    start = time.time()
    windows = get_windows(genome, WINDOW_SIZE, STEP_SIZE, blacklist_path)
    end = time.time()
    print(f"Window creation time: {end - start:.2f} seconds", flush=True)

    # TODO: parallelize this
    for frag_file, sample_name in tqdm(frag_dict.items(), desc = "Processing fragment files") :
        print(f"Processing {sample_name}")
        # Pickle file output
        if PKL_OUTPUT_DIRECTORY is None :
            PKL_OUTPUT_DIRECTORY = f"{OUTPUT_DIRECTORY}/pkl_files"
            os.makedirs(PKL_OUTPUT_DIRECTORY, exist_ok=True)

        if not os.path.exists(PKL_OUTPUT_DIRECTORY):
            os.makedirs(PKL_OUTPUT_DIRECTORY, exist_ok=True)

        # Search for pickle file if already exists
        pickle_out = f"{PKL_OUTPUT_DIRECTORY}/{frag_dict[frag_file]}_windowstats.pkl"
        if Path(pickle_out).exists() :
            print(f"Pickle file found: {pickle_out}")

            with open(pickle_out, "rb") as f:
                data_package = pickle.load(f)
        # Otherwise run pipeline
        else :
            print("Creating cell by window count matrix")
            start = time.time()
            data_package = run_aggregation(frag_file, sample_name, pickle_out, windows, whitelists, N_NEIGHBORS, bgdCN=2, MAKE_TEMP_SAVE=True)
            end = time.time()
            print(f"Cell by window runtime: {end - start:.2f} seconds", flush=True)


        if data_package == None:
            print(f"Count failed for {sample_name}")
            continue
        
        start = time.time()


        # Widen by 1e5
        data_package['wmeta']["Start"] = data_package['wmeta']["Start"] - 100000
        data_package['wmeta']["End"] = data_package['wmeta']["End"] + 100000
        data_package['wmeta']["Start"] = data_package['wmeta']["Start"].clip(lower=0)

        print("Aggregating windows for genes")
        gene_matrix, gene_to_idx = aggregate_genes(genes_pr, data_package)
        col_names = [name for name, idx in sorted(gene_to_idx.items(), key=lambda x: x[1])]
        row_names = [name for name, idx in sorted(data_package['barcodeidx'].items(), key=lambda x: x[1])]

        cnv_df = pd.DataFrame(gene_matrix.T, index=row_names, columns=col_names)

        # Export as TSV
        print("Exporting")
        cnv_df.to_csv(f"{OUTPUT_DIRECTORY}/{frag_dict[frag_file]}_cnv.tsv", sep="\t")

        end = time.time()
        print(f"Gene aggregation runtime {end - start:.2f} seconds", flush=True)
        print("Done")