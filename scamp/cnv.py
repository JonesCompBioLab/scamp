"""
Pure-Python reimplementation of the scATAC CNV pipeline.

This module provides functions to:
- create sliding genome windows from a reference FASTA
- compute GC/AT content for windows
- read 10x-style fragment files
- count fragment insertions per cell per window
- compute background and CNA calls (log2FC, z, pval, padj)
- aggregate window-level copy numbers to gene-level based on a gene BED/GTF

Notes:
- This implementation requires optional dependencies (pyfaidx, pyranges,
  scipy, pandas, numpy). Import errors will provide a friendly message.
- The implementation aims to match the R pipeline numerically but is a
  pragmatic, readable Python-first reimplementation and may differ in
  minor numerical details.

"""

from __future__ import annotations

import gzip
import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# regex to match standard chromosome names
_STANDARD_CHR_RE = re.compile(r'^(chr)?(\d+|X|Y|M|MT)$', re.IGNORECASE)

try:
    import numpy as np
    import pandas as pd
    from pyfaidx import Fasta
    import pysam
    import pyranges as pr
    from scipy import sparse
    from scipy.stats import norm
    from statsmodels.stats.multitest import multipletests
except Exception as exc:  # pragma: no cover - dependency guard
    raise ImportError(
        "scamp.cnv requires optional dependencies (pyfaidx, pysam, pyranges, numpy, pandas, scipy, statsmodels). "
        f"Install scamp[cnv] to use the Python CNV pipeline: {exc}"
    )


def _progress_iter(iterable, description: str | None = None, total: int | None = None, quiet: bool = False):
    """Return an iterator wrapped with a notebook-friendly progress display.

    - If quiet is True, returns the original iterable.
    - If running inside a notebook, tries to use `tqdm.notebook.tqdm`.
    - Otherwise uses `tqdm.tqdm` for console-friendly progress.
    - If `tqdm` is unavailable it returns the original iterable.
    """
    if quiet:
        return iterable

    try:

        from tqdm import tqdm
        return tqdm(iterable, desc=description, total=total)

    except Exception:
        # Either tqdm is not installed or something else went wrong; return raw iterable
        return iterable


def _load_chromosome_sequence(reference_fasta: str, chrom: str) -> str:
    """Load a single chromosome sequence, using pysam if available else pyfaidx.

    Args:
        reference_fasta: Path to reference FASTA.
        chrom: Chromosome name.

    Returns:
        Uppercase sequence string.
    """
    with pysam.FastaFile(reference_fasta) as fasta:
        seq = fasta.fetch(chrom).upper()
    return seq


def _to_ucsc_name(name: str) -> str:
    """Convert chromosome name to UCSC style (e.g., '1' -> 'chr1', 'MT' -> 'chrM').

    If the name does not match standard chromosomes, returns the original name unchanged.
    """
    name_str = str(name)
    m = _STANDARD_CHR_RE.match(name_str)
    if m:
        core = m.group(2)
        if core.upper() in ("M", "MT"):
            return "chrM"
        return f"chr{core}"
    # preserve existing chr-prefixed names
    if name_str.startswith("chr"):
        return name_str
    return name_str


@dataclass
class CNAResult:
    sample_name: str
    windows: pd.DataFrame  # columns: chrom, start, end, name, GC, AT, N
    counts: sparse.csr_matrix  # shape: (n_windows, n_cells)
    cell_names: List[str]
    log2fc: np.ndarray
    padj: np.ndarray


# ------------------ Window creation / GC ------------------ #
def filter_sequences(reference_fasta: Fasta, quiet: bool = False) -> List[Tuple[str, int]]:
    """Filter out non-standard chromosomes and return (name, length) pairs."""
    # Filter out simulated/alternate contigs (keep only standard chromosome names)
    # Accepts names like: '1', '2', ..., 'X', 'Y', 'M', 'MT', optionally prefixed with 'chr'
    standard_re = re.compile(r'^(chr)?(\d+|X|Y|M|MT)$', re.IGNORECASE)
    all_seqs = [(name, len(reference_fasta[name])) for name in reference_fasta.keys()]
    seqs = [s for s in all_seqs if standard_re.match(s[0])]
    if not seqs:
        raise RuntimeError("Reference FASTA contains no standard chromosomes after filtering simulated contigs")
    if len(seqs) < len(all_seqs) and not quiet:
        removed = [n for n, _ in all_seqs if n not in {x[0] for x in seqs}]
        preview = ", ".join(removed[:5])
        print(
            f"Filtered out {len(removed)} non-standard contigs (simulated/alternate). Example: {preview}"
        )
    return seqs

def make_windows(
    reference_fasta: str,
    window_size: int = 3_000_000,
    step_size: int = 1_000_000,
    blacklist_bed: Optional[str] = None,
    min_window_n_fraction: float = 0.999,
    quiet: bool = False,
) -> pd.DataFrame:
    """Create sliding windows across the reference fasta and subtract blacklist.

    This function will:
      - create sliding windows of `window_size` with `step_size`
      - subtract regions in `blacklist_bed` (if provided), producing subwindows
        that retain the original window `name` so they can be grouped later
      - compute GC/AT/N content per resulting subwindow
      - drop subwindows with too much 'N' bases based on
        `min_window_n_fraction` (fraction of non-N required, default 0.999)

    Args:
        quiet: If True, suppress progress bars and other interactive output.

    Returns:
        DataFrame with columns ['Chromosome','Start','End','name','GC','AT','N']
    """
    ref = Fasta(reference_fasta)

    seqs = filter_sequences(ref)
    chroms = [seq[0] for seq in seqs]

    windows = []
    idx = 1
    w_iter = _progress_iter(seqs, description="Making windows", total=len(seqs), quiet=quiet)
    for chrom, length in w_iter:
        # sliding windows

        for start in range(1, length - window_size + 2, step_size):
            end = start + window_size - 1
            name = f"w{idx}"
            windows.append((chrom, start - 1, end, name))  # 0-based start for convenience
            idx += 1

    w_df = pd.DataFrame(windows, columns=["Chromosome", "Start", "End", "name"])

    max_fraction_N = 1 - min_window_n_fraction

    # If a blacklist is provided, subtract it from each window producing pieces.
    if blacklist_bed is not None and Path(blacklist_bed).stat().st_size > 0:
        try:
            blacklist_pr = pr.read_bed(blacklist_bed)
        except Exception as exc:
            raise RuntimeError(f"Failed to read blacklist BED {blacklist_bed}: {exc}")

        windows_pr = pr.PyRanges(w_df[["Chromosome", "Start", "End", "name"]])
        pieces_pr = windows_pr.subtract(blacklist_pr)

        if pieces_pr.empty:
            return pd.DataFrame(columns=["Chromosome", "Start", "End", "name", "GC", "AT", "N"])

        # Keep only required columns; fall through to the unified GC/AT/N computation below
        pieces_df = pieces_pr.df
        pieces_df["name"] = [f"w{i}" for i in range(1, len(pieces_df) + 1)]
        w_df = pieces_df[["Chromosome", "Start", "End", "name"]].reset_index(drop=True)

    # Vectorized GC/AT/N computation: load chromosome sequences on-demand (lazy loading)
    n_rows = len(w_df)
    gcs = np.empty(n_rows, dtype=float)
    ats = np.empty(n_rows, dtype=float)
    ns = np.empty(n_rows, dtype=float)

    # optionally use progress bars (Jupyter-aware)
    outer_iter = _progress_iter(w_df.groupby("Chromosome"), description="Computing GC/AT/N across chromosomes", quiet=quiet)

    for chrom, grp in outer_iter:
        if chrom not in chroms:
            raise RuntimeError(f"Chromosome {chrom} not found in reference FASTA")
        # Load sequence on-demand with optional pysam backend
        seq = _load_chromosome_sequence(reference_fasta, chrom)
        
        starts = grp["Start"].astype(int).to_numpy()
        ends = grp["End"].astype(int).to_numpy()
        idxs = grp.index.to_numpy()

        for idx, s, e in zip(idxs, starts, ends):
            subseq = seq[s:e]
            a = subseq.count("A")
            t = subseq.count("T")
            g = subseq.count("G")
            c = subseq.count("C")
            n = subseq.count("N")
            denom = max(1, a + t + g + c + n)
            gc = (g + c) / denom
            at = (a + t) / denom
            gcs[idx] = gc
            ats[idx] = at
            ns[idx] = 1 - (gc + at)

    w_df["GC"] = gcs
    w_df["AT"] = ats
    w_df["N"] = ns

    # filter out windows with too much N
    keep = w_df.N < max_fraction_N
    w_df = w_df[keep].reset_index(drop=True)

    return w_df


# ------------------ Fragment reading & counting ------------------ #

def read_fragments(fragment_file: str, cell_whitelist: Optional[List[str]] = None, min_frags: int = 100) -> pd.DataFrame:
    """Read a fragment.tsv[.gz] and return DataFrame with columns ['Chromosome','Start','End','cell'].

    Expected format: chrom start end cell barcode (standard 10x fragments.tsv.gz)
    """
    cols = ["Chromosome", "Start", "End", "cell", "reads"]
    compression = "gzip" if str(fragment_file).endswith(".gz") else None
    df = pd.read_csv(
        fragment_file,
        sep="\t",
        header=None,
        names=cols,
        compression=compression,
    )

    # Normalize chromosome names to UCSC style
    df['Chromosome'] = df['Chromosome'].astype(str).apply(_to_ucsc_name)

    if cell_whitelist:
        df = df[df.cell.isin(cell_whitelist)].reset_index(drop=True)

    # keep cells with at least min_frags
    counts = df.cell.value_counts()
    keep_cells = counts[counts >= min_frags].index
    df = df[df.cell.isin(keep_cells)].reset_index(drop=True)
    return df


def count_insertions(
    windows: pd.DataFrame,
    fragments: pd.DataFrame,
    count_ends: bool = True,
    use_fragment_weights: bool = False,
) -> Tuple[sparse.csr_matrix, List[str]]:
    """Count fragment insertions per window per cell (memory-efficient).

    Behavior:
      - If ``count_ends`` is True (default) each fragment contributes two
        point insertions (start and end), matching the R pipeline (counts
        both fragment endpoints).
      - If ``count_ends`` is False, each fragment is treated as a single
        interval and counts are derived by summing a ``reads`` column if
        present or by counting overlaps (one per fragment).

    This implementation processes one chromosome at a time to avoid
    materializing a potentially huge overlap table.

    Returns:
        A sparse matrix of shape (n_windows, n_cells) and ordered list of
        cell names (sorted for determinism).
    """

    # Normalize fragment chromosome names to UCSC to match windows
    fragments = fragments.copy()
    fragments['Chromosome'] = fragments['Chromosome'].astype(str).apply(_to_ucsc_name)

    # Ensure a 'reads' column exists and control whether to treat it as a weight
    if 'reads' not in fragments.columns:
        fragments = fragments.copy()
        fragments['reads'] = 1

    # By default, match R semantics: do not multiply counts by the 'reads' column
    if not use_fragment_weights:
        fragments = fragments.copy()
        fragments['reads'] = 1

    # Precompute window mappings
    window_to_i = {name: i for i, name in enumerate(windows.name)}

    windows_pr = pr.PyRanges(windows[['Chromosome', 'Start', 'End', 'name']])

    # accumulate counts keyed by (window_idx, cell_name)
    counts_acc: dict[tuple[int, str], float] = {}

    # iterate chromosomes present in windows to keep memory bounded
    chroms = sorted(windows['Chromosome'].unique())
    for chrom in chroms:
        # select fragments for this chromosome
        frags_chr = fragments[fragments['Chromosome'] == chrom]
        if frags_chr.shape[0] == 0:
            continue

        # join only windows on this chromosome
        w_chr = windows_pr[windows_pr.Chromosome == chrom]
        if w_chr.empty:
            continue

        if count_ends:
            # create point inserts for both starts and ends
            # represent single-base points as half-open intervals [pos, pos+1)
            starts = frags_chr[['Chromosome', 'Start', 'cell', 'reads']].copy()
            starts['End'] = starts['Start'] + 1
            starts = starts[['Chromosome', 'Start', 'End', 'cell', 'reads']]

            ends = frags_chr[['Chromosome', 'End', 'cell', 'reads']].copy()
            ends['Start'] = ends['End'] - 1
            ends = ends[['Chromosome', 'Start', 'End', 'cell', 'reads']]

            points = pd.concat([starts, ends], ignore_index=True)
            frags_pr = pr.PyRanges(points)
        else:
            frags_pr = pr.PyRanges(frags_chr[['Chromosome', 'Start', 'End', 'cell', 'reads']])

        ov = w_chr.join(frags_pr)
        if ov.empty:
            continue

        ov_df = ov.df
        grouped = ov_df.groupby(['name', 'cell'])['reads'].sum()

        for (wname, cell), val in grouped.items():
            win_i = window_to_i.get(wname)
            if win_i is None:
                # should not happen, but guard against mismatched names
                continue
            counts_acc[(win_i, cell)] = counts_acc.get((win_i, cell), 0.0) + float(val)

    if not counts_acc:
        return sparse.csr_matrix((len(windows), 0)), []

    # build cell index with deterministic sorted order (matches previous behavior)
    unique_cells = sorted({cell for (_, cell) in counts_acc.keys()})
    cell_to_i = {c: i for i, c in enumerate(unique_cells)}

    rows = []
    cols = []
    data = []
    for (win_i, cell), val in counts_acc.items():
        rows.append(win_i)
        cols.append(cell_to_i[cell])
        data.append(val)

    mat = sparse.csr_matrix((np.array(data), (np.array(rows), np.array(cols))),
                             shape=(len(windows), len(unique_cells)), dtype=float)
    return mat, unique_cells


# ------------------ CNA computation ------------------ #

def compute_cna(counts: sparse.csr_matrix, windows: pd.DataFrame, neighbors: int = 100, LFC: float = 1.5, FDR: float = 0.1, force: bool = True) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute background, log2FC, z-score, p-values and FDR-adjusted p-values.

    Returns (cns, log2fc, pvals, padj) shaped (n_windows, n_cells)
    """
    # convert to dense for computations (can be large; consider optimizing later)
    counts_dense = counts.toarray()

    n_windows, n_cells = counts_dense.shape
    bdgMean = np.zeros_like(counts_dense)
    bdgSd = np.zeros_like(counts_dense)
    log2FC = np.zeros_like(counts_dense)
    z = np.zeros_like(counts_dense)
    pval = np.zeros_like(counts_dense)

    gc_vals = windows.GC.to_numpy()

    for i in range(n_windows):
        if i % 100 == 0:
            pass  # could add logging
        # find nearest by GC
        idxNN = np.argsort(np.abs(gc_vals - gc_vals[i]))[: neighbors + 1]
        idxNN = idxNN[idxNN != i]

        if np.any(np.mean(counts_dense[idxNN, :], axis=0) == 0):
            if force:
                # proceed with warning (we don't have a logger here)
                pass
            else:
                raise RuntimeError("Background Mean = 0! Try a higher neighbor count or remove empty cells")

        bdgMean[i, :] = np.mean(counts_dense[idxNN, :], axis=0)
        bdgSd[i, :] = np.std(counts_dense[idxNN, :], axis=0, ddof=1)
        log2FC[i, :] = np.log2( (counts_dense[i, :] + 1e-5) / (bdgMean[i, :] + 1e-5) )

        # z score
        with np.errstate(divide='ignore', invalid='ignore'):
            z[i, :] = (counts_dense[i, :] - bdgMean[i, :]) / bdgSd[i, :]
        pval[i, :] = 2 * norm.cdf(-np.abs(z[i, :]))

    # padj per cell across windows
    padj = np.zeros_like(pval)
    for j in range(n_cells):
        _, adj, _, _ = multipletests(pval[:, j], alpha=FDR, method='fdr_bh')
        padj[:, j] = adj

    cns = 2 * ( 2 ** log2FC ) # convert to copy numbers

    return cns, log2FC, pval, padj


# ------------------ Aggregation to genes ------------------ #

def aggregate_to_genes(counts: np.array, windows: pd.DataFrame, cell_names: List[str], gene_bed: str) -> pd.DataFrame:
    """Aggregate window-level copy numbers to gene bodies.

    Args:
        counts: numpy array (n_windows, n_cells)
        windows: DataFrame with columns Chromosome, Start, End, name
        cell_names: list of cells corresponding to columns of counts
        gene_bed: BED/GTF-like file with columns: chrom, start, end, gene_name

    Returns:
        DataFrame with shape (n_cells, n_genes) (rows=cells, columns=genes)
    """
    # read genes as pyranges
    genes = pr.PyRanges(pd.read_csv(gene_bed, sep='\t', names=['Chromosome', 'Start', 'End', 'strand', 'gene_name']))
    windows_pr = pr.PyRanges(windows[["Chromosome", "Start", "End", "name"]])

    joined = windows_pr.join(genes)
    jdf = joined.df
    if jdf.empty:
        return pd.DataFrame()

    # map each gene to windows
    gene_to_windows = jdf.groupby('gene_name').name.apply(list)  # name_b is gene name column from BED

    # for each gene, average windows across rows
    gene_vals = {}
    for gene, win_names in gene_to_windows.items():
        idxs = [windows.index[windows.name == w].tolist()[0] for w in win_names if w in set(windows.name)]
        if len(idxs) == 0:
            continue
        vals = counts[idxs, :]
        gene_vals[gene] = np.mean(vals, axis=0)

    genes_sorted = sorted(gene_vals.keys())
    df = pd.DataFrame({g: gene_vals[g] for g in genes_sorted}, index=cell_names)
    df.index.name = 'cell'
    return df


# ------------------ High-level pipeline helpers ------------------ #

def run_python_pipeline(
    fragment_directory: str,
    output_directory: str,
    reference_fasta: str,
    gene_bed: str,
    cell_whitelist: Optional[List[str]] = None,
    blacklist_bed: Optional[str] = None,
    window_size: int = 3_000_000,
    step_size: int = 1_000_000,
    n_neighbors: int = 200,
    min_frags: int = 100,
    sample_name: Optional[str] = None,
    quiet: bool = False,
):
    """Run full CNV pipeline on a directory of fragment files.

    Writes gene-level TSV files into `output_directory`.
    Returns list of written gene-level TSV files.
    """
    out_dir = Path(output_directory)
    out_dir.mkdir(parents=True, exist_ok=True)

    ref = reference_fasta

    windows = make_windows(ref, blacklist_bed=blacklist_bed, window_size=window_size, step_size=step_size, quiet=quiet)

    written = []
    frag_files = list(Path(fragment_directory).glob("**/*fragments.tsv.gz"))

    for frag in frag_files:
        sample_name = frag.stem.split('_')[0]

        print(f"Processing {frag.stem} sample with name {sample_name}...")
        if sample_name.lower() in {"atac", "fragments"}:
            warnings.warn(
                f"Detected sample name '{sample_name}' which is one of the reserved names ('atac','fragments'). "
                "This may indicate that fragment files are not named properly and can lead to cell barcodes being misnamed.",
                UserWarning,
            )

        if cell_whitelist:
            cellbc_to_keep = pd.read_csv(cell_whitelist, header=None).values

        frag_df = read_fragments(str(frag), cell_whitelist=cellbc_to_keep, min_frags=min_frags)

        counts, cells = count_insertions(windows, frag_df)
        # keep only windows that have any counts
        if counts.shape[1] == 0:
            continue
        cns, log2fc, pval, padj = compute_cna(counts, windows, neighbors=n_neighbors)
        gene_df = aggregate_to_genes(cns, windows, cells, gene_bed)

        if gene_df.empty:
            continue

        if sample_name:
            gene_df.index = [f'{sample_name}#{c}' for c in gene_df.index]

        out_name = out_dir / f"{sample_name}.all.tsv"
        gene_df.to_csv(out_name, sep='\t', index=False)
        written.append(str(out_name))

    return written
