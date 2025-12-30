import numpy as np
import pandas as pd
import pyranges as pr
from scipy import sparse

from scamp.cnv import count_insertions


def _join_counts_reference(windows: pd.DataFrame, fragments: pd.DataFrame, count_ends: bool = True):
    w_pr = pr.PyRanges(windows)
    if not count_ends:
        f_pr = pr.PyRanges(fragments)
        ov = w_pr.join(f_pr)
        if ov.empty:
            return sparse.csr_matrix((len(windows), 0)), []
        ov_df = ov.df
        grouped = ov_df.groupby(["name", "cell"]).size().reset_index(name="reads")
    else:
        # simulate R behaviour: count both fragment start and end as insertions
        starts = fragments[["Chromosome", "Start", "Start", "cell", 'reads']].copy()
        starts.columns = ["Chromosome", "Start", "End", "cell", 'reads']
        ends = fragments[["Chromosome", "End", "End", "cell", 'reads']].copy()
        ends.columns = ["Chromosome", "Start", "End", "cell", 'reads']
        pts = pd.concat([starts, ends], ignore_index=True)
        f_pr = pr.PyRanges(pts)
        ov = w_pr.join(f_pr)
        if ov.empty:
            return sparse.csr_matrix((len(windows), 0)), []
        ov_df = ov.df
        grouped = ov_df.groupby(["name", "cell"])['reads'].sum().reset_index(name='reads')

    cells = sorted(grouped.cell.unique())
    cell_to_i = {c: i for i, c in enumerate(cells)}
    window_to_i = {name: i for i, name in enumerate(windows.name)}

    rows = grouped.name.map(window_to_i).to_numpy()
    cols = grouped.cell.map(cell_to_i).to_numpy()
    data = grouped.reads.to_numpy()

    mat = sparse.csr_matrix((data, (rows, cols)), shape=(len(windows), len(cells)), dtype=float)
    return mat, cells


def test_count_insertions_matches_join():
    windows = pd.DataFrame([
        {"Chromosome": "chr1", "Start": 0, "End": 100, "name": "w1"},
        {"Chromosome": "chr1", "Start": 100, "End": 200, "name": "w2"},
        {"Chromosome": "chr2", "Start": 0, "End": 150, "name": "w3"},
    ])

    fragments = pd.DataFrame([
        {"Chromosome": "chr1", "Start": 10, "End": 20, "cell": "A", 'reads': 1},
        {"Chromosome": "chr1", "Start": 15, "End": 25, "cell": "A", 'reads': 1},
        {"Chromosome": "chr1", "Start": 110, "End": 115, "cell": "B", 'reads': 1},
        {"Chromosome": "chr2", "Start": 10, "End": 80, "cell": "A", 'reads': 1},
        {"Chromosome": "chr2", "Start": 100, "End": 120, "cell": "C", 'reads': 1},
    ])

    mat_stream, cells_stream = count_insertions(windows, fragments)
    mat_join, cells_join = _join_counts_reference(windows, fragments, count_ends=True)

    assert cells_stream == cells_join
    np.testing.assert_array_equal(mat_stream.toarray(), mat_join.toarray())


def test_count_insertions_with_flag_matches_interval_join():
    # when count_ends=False, count per-fragment (interval) should match join-based reference
    windows = pd.DataFrame([
        {"Chromosome": "chr1", "Start": 0, "End": 100, "name": "w1"},
        {"Chromosome": "chr1", "Start": 100, "End": 200, "name": "w2"},
    ])

    fragments = pd.DataFrame([
        {"Chromosome": "chr1", "Start": 10, "End": 20, "cell": "A", 'reads': 1},
        {"Chromosome": "chr1", "Start": 150, "End": 160, "cell": "B", 'reads': 1},
    ])

    mat_stream, cells_stream = count_insertions(windows, fragments, count_ends=False)
    mat_join, cells_join = _join_counts_reference(windows, fragments, count_ends=False)

    assert cells_stream == cells_join
    np.testing.assert_array_equal(mat_stream.toarray(), mat_join.toarray())


def test_count_insertions_empty_fragments():
    windows = pd.DataFrame([
        {"Chromosome": "chr1", "Start": 0, "End": 100, "name": "w1"},
    ])
    fragments = pd.DataFrame(columns=["Chromosome", "Start", "End", "cell"])

    mat, cells = count_insertions(windows, fragments)
    assert mat.shape == (1, 0)
    assert cells == []
