import os
import numpy as np
from pathlib import Path

import pandas as pd
import pytest
from scipy import sparse


from scamp import cnv


def _write_fasta(fasta_path: Path):
    # chr1 length 1500, chr2 length 1000
    with open(fasta_path, "w") as fh:
        fh.write(">chr1\n")
        fh.write("A" * 300 + "G" * 300 + "C" * 300 + "T" * 600 + "\n")
        fh.write(">chr2\n")
        fh.write("A" * 500 + "G" * 500 + "\n")

def _write_sim_fasta(fp: Path):
    with open(fp, "w") as fh:
        fh.write(">chr1\n")
        fh.write("A" * 1000 + "\n")
        fh.write(">chrUn_random1\n")
        fh.write("A" * 500 + "\n")
        fh.write(">GL000219.1\n")
        fh.write("A" * 200 + "\n")

def test_make_windows_with_blacklist(tmp_path: Path):
    fasta = tmp_path / "ref.fa"
    _write_fasta(fasta)

    # blacklist that overlaps chr1 windows 0-500 and 500-1000
    blacklist = tmp_path / "blacklist.bed"
    with open(blacklist, "w") as fh:
        fh.write("chr1\t200\t800\n")

    # run with window size 500, step 500
    windows = cnv.make_windows(str(fasta), window_size=500, step_size=500, blacklist_bed=str(blacklist))

    # For chr1 we expect three pieces: ~200 (0-200), ~200 (800-1000), and full 500 (1000-1500)
    chr1 = windows[windows.Chromosome == "chr1"].copy()
    assert len(chr1) == 3

    lens = (chr1.End - chr1.Start).to_numpy()
    lens_sorted = np.sort(lens)
    # expected lengths approximately 200, 200, 500
    assert np.isclose(lens_sorted[0], 200, atol=2)
    assert np.isclose(lens_sorted[1], 200, atol=2)
    assert np.isclose(lens_sorted[2], 500, atol=2)

    # Check GC+AT+N sums to ~1 for all returned windows
    sums = (windows["GC"] + windows["AT"] + windows["N"]).to_numpy()
    assert np.allclose(sums, 1.0, atol=1e-6)


def test_make_windows_without_blacklist(tmp_path: Path):
    fasta = tmp_path / "ref2.fa"
    _write_fasta(fasta)

    windows = cnv.make_windows(str(fasta), window_size=500, step_size=500, blacklist_bed=None)

    # chr1 length 1500 -> 3 windows, chr2 length 1000 -> 2 windows => total 5
    assert len(windows) == 5


def test_empty_blacklist_matches_none(tmp_path: Path):
    fasta = tmp_path / "ref3.fa"
    _write_fasta(fasta)

    empty_blacklist = tmp_path / "empty_blacklist.bed"
    empty_blacklist.write_text("")

    w_empty = cnv.make_windows(str(fasta), window_size=500, step_size=500, blacklist_bed=str(empty_blacklist))
    w_none = cnv.make_windows(str(fasta), window_size=500, step_size=500, blacklist_bed=None)

    # Should be identical (same number of windows and the same start/end coordinates)
    assert len(w_empty) == len(w_none)
    assert all(w_empty[["Chromosome","Start","End","name"]].reset_index(drop=True)
               == w_none[["Chromosome","Start","End","name"]].reset_index(drop=True))


def test_blacklist_removes_entire_chromosome(tmp_path: Path):
    fasta = tmp_path / "ref4.fa"
    _write_fasta(fasta)

    blacklist = tmp_path / "blk_full_chr2.bed"
    # Remove chr2 entirely
    with open(blacklist, "w") as fh:
        fh.write("chr2\t0\t1000\n")

    windows = cnv.make_windows(str(fasta), window_size=500, step_size=500, blacklist_bed=str(blacklist))

    # No chr2 windows should remain
    assert len(windows[windows.Chromosome == "chr2"]) == 0
    # chr1 should still have 3 windows
    assert len(windows[windows.Chromosome == "chr1"]) == 3


def test_n_fraction_filtering(tmp_path: Path):
    # Write a fasta with a large N region in chr1
    fasta = tmp_path / "ref_n.fa"
    with open(fasta, "w") as fh:
        fh.write(">chr1\n")
        # 200 A, 400 N, 900 T -> 1500
        fh.write("A" * 200 + "N" * 400 + "T" * 900 + "\n")
        fh.write(">chr2\n")
        fh.write("A" * 500 + "G" * 500 + "\n")

    # With a strict min_window_n_fraction (0.9), windows with >10% N will be dropped
    w_strict = cnv.make_windows(str(fasta), window_size=500, step_size=500, blacklist_bed=None, min_window_n_fraction=0.9)
    # With no restriction, we get more windows
    w_relaxed = cnv.make_windows(str(fasta), window_size=500, step_size=500, blacklist_bed=None, min_window_n_fraction=0.0)

    assert len(w_strict) < len(w_relaxed)
    # Ensure remaining windows have low N fraction
    assert (w_strict["N"] < 0.1).all()


def test_simulated_contigs_filtered(tmp_path: Path):
    fasta = tmp_path / "ref_sim.fa"
    _write_sim_fasta(fasta)

    windows = cnv.make_windows(str(fasta), window_size=500, step_size=500, blacklist_bed=None, quiet=True)

    # All windows should be on chr1 only (simulated contigs filtered out)
    chroms = set(windows.Chromosome.tolist())
    assert chroms == {"chr1"}


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

    mat, cells = cnv.count_insertions(windows, frags)

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



if __name__ == "__main__":
    pytest.main([__file__])
