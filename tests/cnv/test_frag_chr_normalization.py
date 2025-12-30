import pandas as pd
from scipy import sparse

from scamp.cnv import count_insertions


def test_count_insertions_normalizes_chromosome_names():
    # two windows: chr1 and chrM
    windows = pd.DataFrame(
        {
            "Chromosome": ["chr1", "chrM"],
            "Start": [0, 0],
            "End": [1000, 1000],
            "name": ["w1", "w2"],
        }
    )

    # fragments use non-UCSC names: '1' and 'MT'
    frags = pd.DataFrame(
        {
            "Chromosome": ["1", "MT"],
            "Start": [10, 20],
            "End": [11, 21],
            "cell": ["cellA", "cellB"],
            "reads": [1, 2],
        }
    )

    mat, cells = count_insertions(windows, frags)

    # We expect 2 cells and 2 windows
    assert mat.shape == (2, 2)
    assert set(cells) == {"cellA", "cellB"}

    # Check that counts are placed in expected windows
    # convert to dense for assertion
    dense = mat.toarray()
    # Find indices
    cell_to_idx = {c: i for i, c in enumerate(cells)}

    # fragment 0 (cellA) overlaps window 0 (w1)
    assert dense[0, cell_to_idx["cellA"]] >= 1
    # fragment 1 (cellB) overlaps window 1 (w2) with count 2
    assert dense[1, cell_to_idx["cellB"]] >= 2
