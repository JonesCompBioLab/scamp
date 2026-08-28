import gzip

import numpy as np
import pandas as pd
import pyranges as pr
import pytest

from scamp.atac_cnv import cnv_utils


class FakeChromosome:
    def __init__(self, sequence):
        self.sequence = sequence

    def __getitem__(self, key):
        return self.sequence[key]


class FakeGenome:
    def __init__(self, sequences):
        self.sequences = {
            chrom: FakeChromosome(sequence)
            for chrom, sequence in sequences.items()
        }

    def keys(self):
        return self.sequences.keys()

    def __getitem__(self, chrom):
        return self.sequences[chrom]


def _prefix_sums(sequence):
    seq_array = np.frombuffer(sequence.encode(), dtype="S1")
    return {
        "GC": np.cumsum((seq_array == b"G") | (seq_array == b"C")),
        "AT": np.cumsum((seq_array == b"A") | (seq_array == b"T")),
    }


def _write_fragments(path, rows):
    with gzip.open(path, "wt") as handle:
        handle.write("\n".join(rows))
        handle.write("\n")


@pytest.fixture
def simple_windows():
    return pd.DataFrame(
        [
            {
                "Chromosome": "chr1",
                "Start": 0,
                "End": 100,
                "window_id": "chr1:0-100",
                "tile_name": "chr1-0:100",
            }
        ]
    )

def test_subtract_hits_all_overlap_cases():
    main = pr.PyRanges(
        pd.DataFrame(
            [
                ["chr1", 0, 10, "deleted"],
                ["chr1", 20, 30, "trim_start"],
                ["chr1", 40, 50, "trim_end"],
                ["chr1", 60, 80, "split"],
                ["chr1", 90, 100, "untouched"],
            ],
            columns=["Chromosome", "Start", "End", "window_id"],
        )
    )
    blacklist = pr.PyRanges(
        pd.DataFrame(
            [
                ["chr1", 0, 10],
                ["chr1", 15, 25],
                ["chr1", 45, 55],
                ["chr1", 65, 75],
            ],
            columns=["Chromosome", "Start", "End"],
        )
    )

    result = cnv_utils.subtract(main, blacklist).df

    assert set(result.itertuples(index=False, name=None)) == {
        ("chr1", 25, 30, "trim_start"),
        ("chr1", 40, 45, "trim_end"),
        ("chr1", 60, 65, "split"),
        ("chr1", 75, 80, "split"),
        ("chr1", 90, 100, "untouched"),
    }


def test_make_windows_slides_standard_chromosomes_and_applies_blacklist(tmp_path):
    genome = FakeGenome(
        {
            "chr1": "ACGT" * 1500,
            "chr2": "G" * 6000,
            "chrM": "A" * 6000,
        }
    )
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr1\t1000\t2200\n")

    result = cnv_utils.make_windows(
        genome, blacklist, window_size=3000, sliding_size=2000
    )

    segments = set(
        result[["Chromosome", "Start", "End", "window_id"]].itertuples(
            index=False, name=None
        )
    )
    assert segments == {
        ("chr1", 0, 1000, "chr1:0-3000"),
        ("chr1", 2200, 3000, "chr1:0-3000"),
        ("chr1", 2200, 5000, "chr1:2000-5000"),
        ("chr2", 0, 3000, "chr2:0-3000"),
        ("chr2", 2000, 5000, "chr2:2000-5000"),
    }

    assert "chrM" not in set(result["Chromosome"])
    chr1 = result[result["Chromosome"] == "chr1"]
    np.testing.assert_allclose(chr1["gc_count"].to_numpy(), [1400, 400, 500])
    np.testing.assert_allclose(chr1["at_count"].to_numpy(), [1400, 400, 500])

    chr2 = result[result["Chromosome"] == "chr2"]
    np.testing.assert_allclose(chr2["gc_count"].to_numpy(), 3000)
    np.testing.assert_allclose(chr2["at_count"].to_numpy(), 0)


def test_get_whitelists_parses_sample_prefixed_barcodes(tmp_path):
    whitelist = tmp_path / "whitelist.txt"
    whitelist.write_text("sample_a#AAAC\nsample_a#TTTG\nsample_b#CCCC\n")

    result = cnv_utils.get_whitelists(whitelist)

    assert result == {
        "sample_a": ["AAAC", "TTTG"],
        "sample_b": ["CCCC"],
    }


def test_get_whitelists_rejects_unscoped_barcodes(tmp_path):
    whitelist = tmp_path / "bad_whitelist.txt"
    whitelist.write_text("AAAC\n")

    assert cnv_utils.get_whitelists(whitelist) is None


def test_create_cellxwindows_counts_four_column_fragments(tmp_path, simple_windows):
    fragments = tmp_path / "fragments.tsv.gz"
    _write_fragments(
        fragments,
        [
            "chr1\t10\t90\tbc1",
            "chr1\t20\t80\tbc2",
            "chrX\t10\t90\tignored",
        ],
    )
    out_log = []

    result = cnv_utils.create_cellxwindows(
        out_log,
        fragments,
        "sample",
        simple_windows,
        whitelists=None,
        minFrags=1,
        batch_size=1,
    )

    assert set(result.index) == {"bc1", "bc2"}
    assert result.loc["bc1", "chr1-0:100"] == 2
    assert result.loc["bc2", "chr1-0:100"] == 2
    assert "Unique barcodes after filtering: 2" in out_log


def test_create_cellxwindows_accepts_five_column_fragments_and_whitelist(
    tmp_path, simple_windows
):
    fragments = tmp_path / "fragments.tsv.gz"
    _write_fragments(
        fragments,
        [
            "1\t10\t90\tbc1\t3",
            "1\t20\t80\tbc2\t1",
        ],
    )

    result = cnv_utils.create_cellxwindows(
        [],
        fragments,
        "sample",
        simple_windows,
        whitelists={"sample": ["bc1"]},
        minFrags=1,
        batch_size=2,
    )

    assert list(result.index) == ["bc1"]
    assert result.loc["bc1", "chr1-0:100"] == 2


def test_create_cellxwindows_reports_when_no_barcodes_pass_minfrags(
    tmp_path, simple_windows
):
    fragments = tmp_path / "fragments.tsv.gz"
    _write_fragments(fragments, ["chr1\t10\t90\tbc1"])
    out_log = []

    result = cnv_utils.create_cellxwindows(
        out_log,
        fragments,
        "sample",
        simple_windows,
        whitelists=None,
        minFrags=2,
        batch_size=1,
    )

    assert result is None
    assert f"No barcodes passed minFrags (2) in {fragments}" in out_log


def test_create_cellxwindows_rejects_unexpected_fragment_columns(
    tmp_path, simple_windows
):
    fragments = tmp_path / "bad_fragments.tsv.gz"
    _write_fragments(fragments, ["chr1\t10\t90"])

    with pytest.raises(ValueError, match="Unexpected number of columns"):
        cnv_utils.create_cellxwindows(
            [],
            fragments,
            "sample",
            simple_windows,
            whitelists=None,
            minFrags=1,
            batch_size=1,
        )



@pytest.fixture
def fractured_windows():
    return pd.DataFrame(
        [
            {
                "Chromosome": "chr1",
                "Start": 0,
                "End": 30,
                "window_id": "chr1:0-100",
                "tile_name": "chr1:0-30",
                "gc_count" : 10,
                "at_count" : 5,
                "length" : 30
            },
            {
                "Chromosome": "chr1",
                "Start": 35,
                "End": 100,
                "window_id": "chr1:0-100",
                "tile_name": "chr1:80-100",
                "gc_count" : 4,
                "at_count" : 7,
                "length" : 20

            }
        ]
    )



@pytest.fixture
def cellxwindows_df():
    return pd.DataFrame(
        [
            {
                "Cell": "CELL1",
                "chr1:0-30" : 5,
                "chr1:80-100": 4,
            },
            {
                "Cell": "CELL2",
                "chr1:0-30" : 2,
                "chr1:80-100": 12,
            }
        ], index=["CELL1", "CELL2"]
    )

def test_combine_split_windows(fractured_windows, cellxwindows_df) :
    windows_agg, cxw_agg = cnv_utils.combine_split_windows(
        fractured_windows, 100, cellxwindows_df
    )
    assert len(windows_agg) == 1
    assert windows_agg["window_id"].iloc[0] == "chr1:0-100"
    assert windows_agg["gc_count"].iloc[0] == 14
    assert windows_agg["GC_fraction"].iloc[0] == 0.14
    assert windows_agg["N_fraction"].iloc[0] == 0.24
    assert cxw_agg["chr1:0-100"].iloc[0] == 9
    assert cxw_agg["chr1:0-100"].iloc[1] == 14

