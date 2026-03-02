from pathlib import Path
import pandas as pd
import numpy as np
import pyranges as pr
import os
import urllib.request
from pyfaidx import Fasta
from tqdm import tqdm
from collections import defaultdict
import pickle
from scipy.sparse import csr_matrix
import time

'''
Download fasta for hg38 if not already present
Returns Fasta genome
'''
def get_hg38_fasta(cache_dir=f"../reference/.fa_genomes"):
    # hg38 location
    HG38_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"

    cache_dir = os.path.expanduser(cache_dir)
    os.makedirs(cache_dir, exist_ok=True)

    fa_gz = os.path.join(cache_dir, "hg38.fa.gz")
    fa = os.path.join(cache_dir, "hg38.fa")

    # download fasta file as needed
    if not os.path.exists(fa):
        print("Downloading hg38 genome (one-time)...")
        urllib.request.urlretrieve(HG38_URL, fa_gz)

        import gzip, shutil
        with gzip.open(fa_gz, "rb") as f_in, open(fa, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    return Fasta(fa, as_raw=True, sequence_always_upper=True)




'''
Create dict of fragments and file names
Assumes files are formatted as [sample_name]-atac_fragments.tsv.gz
Returns dict: {fragment file : sample name}
'''
def read_frag_files(data_dir, fragment_key) :
    data_dir = Path(data_dir).resolve()
    frag_dict = {}


    # If specified
    if not fragment_key == None :
        with open(fragment_key, "r") as f:
            for line in f :
                frag_dict[f"{data_dir}/{line.strip().split('\t')[0]}"] = frag_dict[line.strip().split('\t')[1]]
        return frag_dict

    # search for .tsv.gz, remove the ending from sample name
    for f in data_dir.glob("*.tsv.gz") :
        sample = f.name.replace("-atac_fragments.tsv.gz", "")
        frag_dict[str(f)] = sample

    return frag_dict


'''
Make pyranges windows, not yet accounting for blacklist
Returns prefix sums for GC and AT content, and windows
'''
def make_windows(genome, REFERENCE_BLACKLIST, window_size = 3000000, sliding_size = 1000000) :
    # Exclude experimental chromosomes
    standard_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    chroms = [c for c in genome.keys() if c in standard_chroms]

    # Set up blacklist
    blacklist = pr.read_bed(f"{REFERENCE_BLACKLIST}")
    # Remove small sections of blacklist
    blacklist = blacklist[ (blacklist.End - blacklist.Start) > 1000 ]

    # {chr:{'GC':[...], 'AT':[...]}}
    # N content is calculated from the other two
    # prefix_sums = {}

    # list of window dicts
    final_dfs = []

    print("Computing windows...")
    for chr in chroms :
        print(f"{chr}")
        seq = genome[chr][:]
        
        # Faster than a string list
        seq_array = np.frombuffer(seq.encode(), dtype='S1')

        gc_flags = (seq_array == b'G') | (seq_array == b'C')
        at_flags = (seq_array == b'A') | (seq_array == b'T')

        prefix_sums = {'GC' : np.cumsum(gc_flags), 'AT' : np.cumsum(at_flags)}
        chrom_len = len(prefix_sums["GC"])

        rows = []
        for start in range(0, chrom_len - window_size + 1, sliding_size):
            end = start + window_size
            rows.append({
                "Chromosome": chr,
                "Start": start,
                "End": end,
                "window_id": f"{chr}:{start}-{end}"
            })

        temp_windows = pr.PyRanges(pd.DataFrame(rows))

        # Subtract blacklist
        temp_windows_m_blacklist = subtract(temp_windows, blacklist)
        temp_windows_m_blacklist = temp_windows_m_blacklist.df.copy()

        # Calculate GC, AT fractions
        temp_windows_m_blacklist_fracs = calculate_fractions(temp_windows_m_blacklist, prefix_sums)
        final_dfs.append(temp_windows_m_blacklist_fracs)

    windows_m_blacklist_fracs = pd.concat(final_dfs, ignore_index=True)

    return windows_m_blacklist_fracs

'''
Read and extract whitelist from whitelist file
Whitelist file should be [sample_name]#barcode
Returns: a dict, {file_name : list_of_barcodes}
'''
def get_whitelists(WHITELIST_FILE) :
    if WHITELIST_FILE :
        whitelist_full = []
        with open(WHITELIST_FILE, "r") as f:
            whitelist_full = [line.strip() for line in f]

        whitelists = defaultdict(list)
        # sample_name : barcodes
        for line in whitelist_full :

            # Ensure format
            if not '#' in line :
                print("Warning: whitelist format should be [sample_name]#[barcode], using all cells as whitelist")
                return None
            whitelists[line.split('#')[0]].append(line.split('#')[1])
            
    else :
        whitelists = None

    return whitelists

'''
Subtracts two pyranges
main: pr to subtract from
to_subtract: what to subtract
Returns pyrange subtracted objected
'''
def subtract(main, to_subtract) :
    # Sorted start/end lists per chromosome
    intervals = defaultdict(list)

    for chrom, start, end, wid in zip(
        main.df["Chromosome"],
        main.df["Start"],
        main.df["End"],
        main.df["window_id"],
    ):
        intervals[chrom].append((start, end, wid))

    for row_num, row in to_subtract.df.iterrows() :
        start, end, chr = row["Start"], row["End"], row["Chromosome"]

        # Go from the back to remove from list
        for idx in range(len(intervals[chr]) - 1, -1, -1) :
            # [start, end)
            # Completely inside blacklist, delete
            if start <= intervals[chr][idx][0] and end >= intervals[chr][idx][1] :
                del intervals[chr][idx]
            # Covers start (can ignore the completely inside case)
            elif start <= intervals[chr][idx][0] and end > intervals[chr][idx][0] :
                intervals[chr].append((end, intervals[chr][idx][1], intervals[chr][idx][2]))
                del intervals[chr][idx]
            # Covers end
            elif start < intervals[chr][idx][1] and end >= intervals[chr][idx][1] :
                intervals[chr].append((intervals[chr][idx][0], start, intervals[chr][idx][2]))
                del intervals[chr][idx]
            # blacklist inside int (split into two)
            elif start < intervals[chr][idx][1] and end > intervals[chr][idx][0] :
                # add two new windows
                intervals[chr].append((end, intervals[chr][idx][1], intervals[chr][idx][2])) 
                intervals[chr].append((intervals[chr][idx][0], start, intervals[chr][idx][2]))
                del intervals[chr][idx]

    # Turn back into pyranges
    rows = []

    for chrom, intervals in intervals.items():
        for start, end, window_id in intervals:
            rows.append([chrom, start, end, window_id])

    df = pd.DataFrame(rows, columns=["Chromosome", "Start", "End", "window_id"])
    gr = pr.PyRanges(df)

    return gr

'''
Get GC, AT fraction annotations for windows
Returns pyranges for updated windows
'''
def calculate_fractions(windows, prefix_sums):
    gc_fractions = []
    at_fractions = []
    n_fractions = []
    for i, window in windows.iterrows() :

        # Calculate the GC, AT, and N counts in the window using prefix sums
        if window['Start'] > 0 :
            gc_count = prefix_sums['GC'][window['End']-1] - prefix_sums['GC'][window['Start']-1]
            at_count = prefix_sums['AT'][window['End']-1] - prefix_sums['AT'][window['Start']-1]
        else :
            gc_count = prefix_sums['GC'][window['End']-1]
            at_count = prefix_sums['AT'][window['End']-1]

        total_bases = (window['End'] - window['Start']) 

        gc_fractions.append(gc_count/total_bases)
        at_fractions.append(at_count/total_bases)
        n_fractions.append(1 - gc_count/total_bases - at_count/total_bases)
    
    windows['GC_fraction'] = gc_fractions
    windows['AT_fraction'] = at_fractions
    windows['N_fraction'] = n_fractions

    return windows

def weighted_avg(subdf, col):
    # weighted by width
    return (subdf[col] * subdf['width']).sum() / subdf['width'].sum()


'''
Creates window file if not yet made, otherwise opens the windows file and cleans it
The name is created based on window size and step size, so inputting different values will make a new file
Returns annotated windows
'''
def get_windows(genome, WINDOW_SIZE, STEP_SIZE, REFERENCE_BLACKLIST) :
    # Create windows file if not yet created
    if not os.path.exists(f"{REFERENCE_BLACKLIST}_{WINDOW_SIZE}_{STEP_SIZE}.tsv") :
        print("Windows not found, creating")
        start = time.time()

        # Calculate prefix sums
        windows_m_blacklist_fracs = make_windows(genome, REFERENCE_BLACKLIST, WINDOW_SIZE, STEP_SIZE)

        # Export
        windows_m_blacklist_fracs.to_csv(f'{REFERENCE_BLACKLIST}_{WINDOW_SIZE}_{STEP_SIZE}.tsv', sep = '\t', index = False)
        windows_m_blacklist_fracs = pd.read_csv(f"{REFERENCE_BLACKLIST}_{WINDOW_SIZE}_{STEP_SIZE}.tsv", sep = '\t')
        end = time.time()
        print(f"Window creation time: {end - start:.2f} seconds", flush=True)

    # If file already created, just read it in
    else :
        windows_m_blacklist_fracs = pd.read_csv(f"{REFERENCE_BLACKLIST}_{WINDOW_SIZE}_{STEP_SIZE}.tsv", sep = '\t')

    # Recombine windows
    windows_m_blacklist_fracs['width'] = windows_m_blacklist_fracs['End'] - windows_m_blacklist_fracs['Start']
    windows_m_blacklist_fracs = windows_m_blacklist_fracs.groupby('window_id').apply(
    lambda x: pd.Series({
        'Chromosome': x['Chromosome'].iloc[0],
        'Start': x['Start'].min(),
        'End': x['End'].max(),
        'GC_fraction': weighted_avg(x, 'GC_fraction'),
        'AT_fraction': weighted_avg(x, 'AT_fraction'),
        'N_fraction': weighted_avg(x, 'N_fraction')
        })
    ).reset_index()

    # Remove windows with large fraction as N
    windows_m_blacklist_fracs = windows_m_blacklist_fracs.loc[windows_m_blacklist_fracs["N_fraction"] < 0.001]
    chrs_to_use = [f"chr{i}" for i in range(1, 23)]
    windows_m_blacklist_fracs = windows_m_blacklist_fracs[windows_m_blacklist_fracs["Chromosome"].isin(chrs_to_use)]

    # add tile id
    windows_m_blacklist_fracs["tile_name"] = windows_m_blacklist_fracs["Chromosome"].astype(str) + "-" + windows_m_blacklist_fracs["Start"].astype(str) + ":" + windows_m_blacklist_fracs["End"].astype(str)

    # Drop duplicate rows
    windows_m_blacklist_fracs = windows_m_blacklist_fracs.drop_duplicates(subset="tile_name", keep="first").reset_index(drop=True)

    return windows_m_blacklist_fracs


'''
Creates cell by window matrix
Batch size is the number of points it looks at at a time (twice the number of fragments)
Returns pandas dataframe representation
'''
def create_cellxwindows(out_log, frag_file, sample_name, windows, whitelists, minFrags = 100, batch_size = 1000000) :
    frag_df = pd.read_csv(
        frag_file,
        sep="\t",
        header=None,
        compression="gzip",
        comment="#",
    )

    out_log.append("Frag file loaded")
    
    # Resolve fragment file differences
    ncol = frag_df.shape[1]
    if ncol == 4:
        frag_df.columns = ["Chromosome", "Start", "End", "Barcode"]
    elif ncol == 5:
        frag_df.columns = ["Chromosome", "Start", "End", "Barcode", "Count"]
    else:
        raise ValueError(f"Unexpected number of columns in fragment file {frag_file}: {ncol}")
    
    if whitelists is not None and sample_name in whitelists :
        mask = frag_df["Barcode"].isin(whitelists[sample_name])
        frag_df = frag_df.loc[mask] 
    

    # Remove chromosomes
    chrs_to_use = [f"chr{i}" for i in range(1, 23)]
    frag_df = frag_df[
        frag_df["Chromosome"].isin(chrs_to_use)
    ]

    if len(frag_df) == 0 :
        out_log.append(f"WARNING: no cells passed minFrag & whitelist in {frag_file}, excluding")
        return None

    # Remove minfrags
    out_log.append("Removing based on minfrags")
    barcodes = frag_df['Barcode'].unique()
    barcode_counts = {
        barcode: idx.to_numpy()
        for barcode, idx in frag_df.groupby("Barcode").groups.items()
        if len(idx) >= minFrags
    }
    barcodes = barcode_counts.keys()
    mask = frag_df["Barcode"].isin(barcodes)
    frag_df = frag_df.loc[mask]  
    out_log.append(f"Unique barcodes after processing: {len(barcodes)}")  
    out_log.append(f"Number of fragments: {len(frag_df)}")
    out_log.append(f"Number of windows: {len(windows)}") 

    windows_pr = pr.PyRanges(windows)

    # Fragment starts
    starts_df = frag_df[["Chromosome", "Start", "Barcode"]].copy()
    starts_df["End"] = starts_df["Start"]  # make it a single-point interval
    starts_df["Count"] = 1  # optional

    # Fragment ends
    ends_df = frag_df[["Chromosome", "End", "Barcode"]].copy()
    ends_df = ends_df.rename(columns={"End":"Start"})
    ends_df["End"] = ends_df["Start"]  # single-point interval
    ends_df["Count"] = 1

    tiles = windows['tile_name'].unique()
    cellxwindow_df = np.zeros((len(barcodes), len(tiles)), dtype=int)
    barcode_idx = {b: i for i, b in enumerate(barcodes)}
    tile_idx = {t: i for i, t in enumerate(tiles)}

    # Combine in chunks
    points_df = pd.concat([starts_df, ends_df], axis=0)
    out_log.append("Overlapping...")
    for start in range(0, len(points_df), batch_size):
        # out_log.append(f"Fragments {start/2} to {(start + batch_size)/2}")
        frag_chunk = points_df.iloc[start:start+batch_size]
        points_pr = pr.PyRanges(frag_chunk)
        overlaps = points_pr.join(windows_pr)
        # For each overlap, increment the count
        r = overlaps.df["Barcode"].map(barcode_idx).to_numpy()
        c = overlaps.df["tile_name"].map(tile_idx).to_numpy()

        np.add.at(cellxwindow_df, (r, c), 1)
    
        del overlaps

    # Output
    cellxwindow_df = pd.DataFrame(
        cellxwindow_df,
        index=barcodes,
        columns=tiles
    )
    return cellxwindow_df


'''
Turn gene annotations into right format for pyranges join
Returns correct pyranges
'''
def get_gene_pr(genes_path) :
    genes_df = pd.read_csv(genes_path, sep = '\t')
    genes_df = genes_df.rename(columns={"seqnames" : "Chromosome", "start" : "Start", "end" : "End"})
    return pr.PyRanges(genes_df)


'''
Get counts per cell for each gene
Returns
'''
def aggregate_genes(genes_pr, data_package) :
    windows_pr = pr.PyRanges(data_package["wmeta"])

    overlaps = genes_pr.join(windows_pr)

    overlap_nodup = overlaps.df[["symbol", "window_id"]].drop_duplicates()
    overlap_nodup["row_id"] = overlap_nodup["window_id"].map(data_package["windowidx"]).astype(int)

    genes = overlap_nodup["symbol"].unique()
    gene_to_idx = {g: i for i, g in enumerate(genes)}

    rows = []
    cols = []
    data = []

    # group windows per gene
    for gene, sub in overlap_nodup.groupby("symbol"):
        g_idx = gene_to_idx[gene]
        wrows = sub["row_id"].values

        weight = 1.0 / len(wrows)

        rows.extend([g_idx] * len(wrows))
        cols.extend(wrows)
        data.extend([weight] * len(wrows))

    # Use matrix multiplication to compute
    W = csr_matrix(
        (data, (rows, cols)),
        shape=(len(genes), data_package["CNs"].shape[0])
    )
    gene_matrix = W @ data_package["CNs"] 
    
    return gene_matrix, gene_to_idx


'''
Cell count and aggregation pipeline
'''
def run_aggregation(out_log, frag_file, sample_name, pickle_out, windows, whitelists, neighbors = 200, bgdCN = 2, MAKE_TEMP_SAVE = True) :
    # cell by windows matrix
    cellxwindows_df = create_cellxwindows(out_log, frag_file, sample_name, windows, whitelists)

    if cellxwindows_df is None :
        return None

    # Remove windows where all are empty
    nonzero_windows = cellxwindows_df.columns
    windows_r = windows[
        windows["tile_name"].isin(nonzero_windows)
    ].copy()

    # Collect windows metadata
    # TODO: I don't think the script actually uses the effective length, though the blacklist used is very short so it wouldn't matter
    wmeta = (
        windows_r
        .assign(width=lambda d: d.End - d.Start)
        .groupby("window_id")
        .apply(lambda d: pd.Series({

            "effectiveLength": d["width"].sum(),
            "percentEffectiveLength": 100 * d["width"].sum() / (d.End.max() - d.Start.min()),
            "GC": (d.GC_fraction * d.width).sum() / d.width.sum(),
            "AT": (d.AT_fraction * d.width).sum() / d.width.sum(),
            "N":  (d.N_fraction  * d.width).sum() / d.width.sum(),
            "Chromosome": d.Chromosome.iloc[0],
            "Start": d.Start.min(),
            "End": d.End.max(),
        })).reset_index() 
    )

    # get counts by window ids
    tile2window = {}
    for i, row in windows_r.iterrows() :
        tile2window[row["tile_name"]] = row["window_id"]

    window_ids = cellxwindows_df.columns.map(tile2window)
    countSummary = cellxwindows_df.T.groupby(window_ids).sum()

    window_ids_final = {id : i for i, id in enumerate(countSummary.index.to_list())}
    cell_barcodes = {barcode : i for i, barcode in enumerate(countSummary.columns.to_list())}
    X = countSummary.to_numpy() 

    gc = wmeta["GC"].to_numpy()

    # Calculate statistics
    bdgMean = np.zeros_like(X, dtype=float)
    log2FC  = np.zeros_like(X, dtype=float)

    for i in range(len(gc)):
        
        # find nearest neighbors by GC
        distances = np.abs(gc[i] - gc)
        nn_idx = np.argsort(distances)[1 : neighbors+1]

        # background
        bg = X[nn_idx, :]
        bg_mean = bg.mean(axis=0)
        bdgMean[i, :] = bg_mean

        # Stats from bg
        log2FC[i, :] = np.log2((X[i, :] + 1e-5) / (bg_mean + 1e-5))

    CNs = bgdCN * (2 ** log2FC)

    # Resulting data
    data_package = {"wmeta" : wmeta, "windows" : windows_r, "windowidx" : window_ids_final, "barcodeidx": cell_barcodes, "CNs" : CNs, "bdgMean" : bdgMean, "counts" : X, "log2FC" : log2FC}
    
    # Export (temporary save)
    if MAKE_TEMP_SAVE :
        with open(pickle_out, "wb") as out :
            pickle.dump(data_package, out)

    return data_package
