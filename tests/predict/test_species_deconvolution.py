from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.cluster.hierarchy import fcluster, linkage

from scamp.predict import species_deconvolution


def _species_cells(n_cells=300, seed=0, overlap=False, zero_cells=0):
    """Copy-numbers of genes on two ecDNA species, A and B.

    Each cell carries `a` copies of species A and `b` of species B. A1, A2, B1 and B2
    are on their species once, A3 and B3 twice, so their copy-number is 2 + a, or
    2 + 2a. With `overlap`, genes AB1-AB3 are on both species (2 + a + b). The first
    `zero_cells` cells have no ecDNA (copy-number 2 everywhere).

    Returns:
        A tuple `(cellbygene_df, a, b)`.
    """
    rng = np.random.default_rng(seed)
    a = rng.gamma(2.0, 5.0, n_cells)
    b = rng.gamma(2.0, 5.0, n_cells)
    data = {}
    for i, copies in enumerate([1, 1, 2]):
        data[f"A{i + 1}"] = 2 + copies * a + rng.normal(0, 0.05, n_cells)
        data[f"B{i + 1}"] = 2 + copies * b + rng.normal(0, 0.05, n_cells)
    if overlap:
        for i in range(3):
            data[f"AB{i + 1}"] = 2 + a + b + rng.normal(0, 0.05, n_cells)
    cellbygene_df = pd.DataFrame(data, index=[f"cell{i}" for i in range(n_cells)]).clip(lower=2)
    if zero_cells:
        cellbygene_df.iloc[:zero_cells] = 2.0
        a[:zero_cells] = 0
        b[:zero_cells] = 0
    return cellbygene_df, a, b


def _groups(species_to_gene):
    """The species' gene sets, ignoring species names and gene order."""
    return {frozenset(genes) for genes in species_to_gene.values()}


def _species_of(species_to_gene, gene):
    return next(species for species, genes in species_to_gene.items() if gene in genes)


TWO_SPECIES = {frozenset({"A1", "A2", "A3"}), frozenset({"B1", "B2", "B3"})}

# Each gene has noise with standard deviation 0.05, so an amount estimated from it is
# within about 6 standard deviations of the truth in every one of the 300 cells
AMOUNT_TOL = 0.2


# remove_duplicates


def test_remove_duplicates_keeps_first_gene_of_each_group():
    cellbygene = pd.DataFrame(
        {"G1": [1.0, 2.0, 3.0], "G2": [5.0, 1.0, 0.0], "G3": [1.0, 2.0, 3.0]}
    )

    df_unique, column_groups = species_deconvolution.remove_duplicates(cellbygene)

    assert list(df_unique.columns) == ["G1", "G2"]
    assert column_groups == {"G1": ["G1", "G3"], "G2": ["G2"]}


def test_remove_duplicates_uses_absolute_tolerance():
    cellbygene = pd.DataFrame(
        {"G1": [1.0, 2.0], "G2": [1.0 + 1e-7, 2.0], "G3": [1.0 + 1e-3, 2.0]}
    )

    _, column_groups = species_deconvolution.remove_duplicates(cellbygene, tol=1e-5)

    assert column_groups == {"G1": ["G1", "G2"], "G3": ["G3"]}


# residualize and fwl_beta


def test_residualize_leaves_residuals_orthogonal_to_the_covariates():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(100, 2))
    Y = rng.normal(size=(100, 3))

    residuals = species_deconvolution.residualize(Y, X)

    assert residuals.shape == Y.shape
    np.testing.assert_allclose(X.T @ residuals, 0, atol=1e-8)


def test_residualize_removes_targets_in_the_span_of_the_covariates():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(50, 2))

    residuals = species_deconvolution.residualize(X @ np.array([3.0, -1.0]), X)

    assert residuals.shape == (50,)
    np.testing.assert_allclose(residuals, 0, atol=1e-8)


def test_fwl_beta_matches_the_full_regression_coefficient():
    rng = np.random.default_rng(0)
    x1, x2 = rng.normal(size=(2, 200))
    y = 2.0 * x1 + 3.0 * x2 + rng.normal(0, 0.1, 200)
    full_coefficients, *_ = np.linalg.lstsq(np.column_stack([x1, x2]), y, rcond=None)

    other = x1.reshape(-1, 1)
    beta = species_deconvolution.fwl_beta(
        species_deconvolution.residualize(x2, other), species_deconvolution.residualize(y, other)
    )

    assert beta == pytest.approx(full_coefficients[1])


# _silhouette_with_extra and _hier_get_k


def test_silhouette_with_extra_scores_a_single_cluster():
    embed = np.corrcoef(_species_cells()[0], rowvar=False)

    silhouette = species_deconvolution._silhouette_with_extra(
        embed, np.ones(embed.shape[0], dtype=int), ddist=1.3, metric="euclidean"
    )

    assert -1 <= silhouette <= 1


def test_hier_get_k_finds_two_species():
    embed = np.corrcoef(_species_cells()[0], rowvar=False)

    assert species_deconvolution._hier_get_k(embed, None, 1.3) == 2


def test_hier_get_k_respects_max_species():
    embed = np.corrcoef(_species_cells()[0], rowvar=False)

    assert species_deconvolution._hier_get_k(embed, 1, 1.3) == 1


def test_hier_get_k_leeway_favours_more_species_and_defaults_to_genes_minus_one():
    # A leeway above the silhouette range makes every number of species a tie, so the
    # largest one checked is chosen: the number of genes (6) minus 1 by default
    embed = np.corrcoef(_species_cells()[0], rowvar=False)

    assert species_deconvolution._hier_get_k(embed, None, 1.3, leeway=2) == 5


# _check_overlap
def _clusters(cellbygene_df, k):
    embed = np.corrcoef(cellbygene_df, rowvar=False)
    return fcluster(linkage(embed, method="average", metric="correlation"), t=k, criterion="maxclust")


def test_check_overlap_needs_at_least_three_clusters():
    cellbygene_df, _, _ = _species_cells(overlap=True)

    assert not species_deconvolution._check_overlap(
        cellbygene_df, _clusters(cellbygene_df, 2), overlap_thresh=0.55
    )


def test_check_overlap_detects_genes_shared_between_species():
    cellbygene_df, _, _ = _species_cells(overlap=True)

    assert species_deconvolution._check_overlap(
        cellbygene_df, _clusters(cellbygene_df, 3), overlap_thresh=0.55
    )


def test_check_overlap_is_false_for_separate_species():
    cellbygene_df, _, _ = _species_cells()

    assert not species_deconvolution._check_overlap(
        cellbygene_df, _clusters(cellbygene_df, 3), overlap_thresh=0.55
    )


# hierarchical_deconvolution
def test_hierarchical_deconvolution_recovers_two_species():
    cellbygene_df, _, _ = _species_cells()

    species_to_gene, cell_by_ecDNA = species_deconvolution.hierarchical_deconvolution(
        cellbygene_df
    )

    assert _groups(species_to_gene) == TWO_SPECIES
    assert set(cell_by_ecDNA.columns) == set(species_to_gene)
    assert list(cell_by_ecDNA.index) == list(cellbygene_df.index)


def test_hierarchical_deconvolution_uses_given_number_of_species():
    cellbygene_df, _, _ = _species_cells()

    species_to_gene, _ = species_deconvolution.hierarchical_deconvolution(
        cellbygene_df, num_ecDNA=1
    )

    assert _groups(species_to_gene) == {frozenset(cellbygene_df.columns)}


def test_hierarchical_deconvolution_amount_averages_genes_near_the_lowest():
    # A3 has twice A1's total (more than 1.3 times), so only A1 and A2 are averaged
    cellbygene_df, a, b = _species_cells()

    species_to_gene, cell_by_ecDNA = species_deconvolution.hierarchical_deconvolution(
        cellbygene_df
    )

    np.testing.assert_allclose(cell_by_ecDNA[_species_of(species_to_gene, "A1")], a, atol=AMOUNT_TOL)
    np.testing.assert_allclose(cell_by_ecDNA[_species_of(species_to_gene, "B1")], b, atol=AMOUNT_TOL)


# naive_deconvolution
def test_naive_deconvolution_recovers_two_species():
    cellbygene_df, a, _ = _species_cells()

    species_to_gene, cell_by_ecDNA = species_deconvolution.naive_deconvolution(cellbygene_df)

    assert _groups(species_to_gene) == TWO_SPECIES
    np.testing.assert_allclose(cell_by_ecDNA[_species_of(species_to_gene, "A1")], a, atol=AMOUNT_TOL)


def test_naive_deconvolution_uses_given_number_of_species():
    cellbygene_df, _, _ = _species_cells()

    species_to_gene, _ = species_deconvolution.naive_deconvolution(cellbygene_df, num_ecDNA=1)

    assert _groups(species_to_gene) == {frozenset(cellbygene_df.columns)}


# _hier_nnls
def test_hier_nnls_recovers_species_and_amounts():
    # A3 counts as two copies of species A, so the NNLS amount is still `a`
    cellbygene_df, a, b = _species_cells()

    species_to_gene, cell_by_ecDNA = species_deconvolution._hier_nnls(cellbygene_df)

    assert _groups(species_to_gene) == TWO_SPECIES
    np.testing.assert_allclose(cell_by_ecDNA[_species_of(species_to_gene, "A1")], a, atol=AMOUNT_TOL)
    np.testing.assert_allclose(cell_by_ecDNA[_species_of(species_to_gene, "B1")], b, atol=AMOUNT_TOL)


# _cNMF_usage
def test_cnmf_usage_rescales_masks_and_names_columns_after_spectra():
    genes = ["G1", "G2", "G3", "G4"]
    species = ["pred_ecDNA_1", "pred_ecDNA_2"]
    # G3 has a little weight in species 1, but is not one of its genes, so it is masked
    spectra_tpm = pd.DataFrame(
        [[2.0, 0.0], [2.0, 0.0], [0.5, 3.0], [0.0, 3.0]], index=genes, columns=species
    )
    species_to_gene = {"pred_ecDNA_1": ["G1", "G2"], "pred_ecDNA_2": ["G3", "G4"]}
    true_usage = np.array([[1.0, 0.0], [2.0, 5.0], [0.0, 4.0]])
    cellbygene = pd.DataFrame(
        true_usage @ np.array([[1.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 1.0]]),
        index=["c1", "c2", "c3"],
        columns=genes,
    )
    usage_df = pd.DataFrame(0.0, index=cellbygene.index, columns=[1, 2])

    result = species_deconvolution._cNMF_usage(
        species, species_to_gene, spectra_tpm, cellbygene, usage_df
    )

    assert list(result.columns) == species
    np.testing.assert_allclose(result.to_numpy(), true_usage, atol=1e-8)


# cNMF_deconvolution
def test_cnmf_deconvolution_recovers_two_species(tmp_path):
    cellbygene_df, a, _ = _species_cells()

    species_to_gene, cell_by_ecDNA = species_deconvolution.cNMF_deconvolution(
        cellbygene_df, "sample", [], n_iter=5, num_ecDNA=2, log_dir=str(tmp_path), seed=1
    )

    assert _groups(species_to_gene) == TWO_SPECIES
    assert list(cell_by_ecDNA.columns) == list(species_to_gene)
    assert list(cell_by_ecDNA.index) == list(cellbygene_df.index)
    assert np.corrcoef(cell_by_ecDNA[_species_of(species_to_gene, "A1")], a)[0, 1] > 0.99


def test_cnmf_deconvolution_chooses_the_number_of_species(tmp_path):
    cellbygene_df, _, _ = _species_cells()

    species_to_gene, _ = species_deconvolution.cNMF_deconvolution(
        cellbygene_df, "sample", [], n_iter=5, max_species=3, log_dir=str(tmp_path), seed=1
    )

    assert _groups(species_to_gene) == TWO_SPECIES


def test_cnmf_deconvolution_names_species_from_one(tmp_path):
    cellbygene_df, _, _ = _species_cells()

    species_to_gene, cell_by_ecDNA = species_deconvolution.cNMF_deconvolution(
        cellbygene_df, "sample", [], n_iter=5, num_ecDNA=3, log_dir=str(tmp_path), seed=1
    )

    expected = [f"pred_ecDNA_{i + 1}" for i in range(len(species_to_gene))]
    assert list(species_to_gene) == expected
    assert list(cell_by_ecDNA.columns) == expected


def test_cnmf_deconvolution_gives_cells_without_ecdna_zero_usage(tmp_path):
    cellbygene_df, _, _ = _species_cells(zero_cells=5)
    out_log = []

    _, cell_by_ecDNA = species_deconvolution.cNMF_deconvolution(
        cellbygene_df, "sample", out_log, n_iter=5, num_ecDNA=2, log_dir=str(tmp_path), seed=1
    )

    assert len(cell_by_ecDNA) == len(cellbygene_df)
    np.testing.assert_array_equal(cell_by_ecDNA.iloc[:5].to_numpy(), 0)


def test_cnmf_deconvolution_removes_its_temporary_folder(tmp_path):
    cellbygene_df, _, _ = _species_cells()

    species_deconvolution.cNMF_deconvolution(
        cellbygene_df, "sample", [], n_iter=5, num_ecDNA=2, log_dir=str(tmp_path), seed=1
    )

    assert not (tmp_path / "sample" / "cnmf_tmp").exists()
    assert (tmp_path / "sample").exists()


def test_cnmf_deconvolution_rejects_more_species_than_max_species(tmp_path):
    cellbygene_df, _, _ = _species_cells()

    with pytest.raises(ValueError, match="larger than max_species"):
        species_deconvolution.cNMF_deconvolution(
            cellbygene_df, "sample", [], num_ecDNA=3, max_species=2, log_dir=str(tmp_path)
        )

    assert not any(Path(tmp_path).iterdir())


# hybrid_deconvolution
def test_hybrid_deconvolution_uses_hierarchical_clusters_without_overlap(tmp_path):
    cellbygene_df, a, _ = _species_cells()
    out_log = []
    log_dir = tmp_path / "cnmf"

    species_to_gene, cell_by_ecDNA = species_deconvolution.hybrid_deconvolution(
        cellbygene_df, out_log, log_dir=str(log_dir), verbose=True
    )

    assert _groups(species_to_gene) == TWO_SPECIES
    np.testing.assert_allclose(cell_by_ecDNA[_species_of(species_to_gene, "A1")], a, atol=AMOUNT_TOL)
    assert "Detected Overlap" not in out_log
    assert not log_dir.exists()


def test_hybrid_deconvolution_runs_cnmf_up_to_the_hierarchical_number_on_overlap(monkeypatch):
    cellbygene_df, _, _ = _species_cells(overlap=True)
    calls = []
    monkeypatch.setattr(
        species_deconvolution,
        "cNMF_deconvolution",
        lambda *args, **kwargs: calls.append(kwargs) or ({}, pd.DataFrame()),
    )
    out_log = []

    species_deconvolution.hybrid_deconvolution(cellbygene_df, out_log, verbose=True)

    assert "Detected Overlap" in out_log
    assert len(calls) == 1
    assert calls[0]["max_species"] == 3
    assert calls[0]["num_ecDNA"] is None


def test_hybrid_deconvolution_passes_given_number_of_species_to_cnmf(monkeypatch):
    cellbygene_df, _, _ = _species_cells(overlap=True)
    calls = []
    monkeypatch.setattr(
        species_deconvolution,
        "cNMF_deconvolution",
        lambda *args, **kwargs: calls.append(kwargs) or ({}, pd.DataFrame()),
    )

    species_deconvolution.hybrid_deconvolution(cellbygene_df, [], num_ecDNA=3)

    assert calls[0]["max_species"] == 3
    assert calls[0]["num_ecDNA"] == 3


def test_hybrid_deconvolution_assigns_every_gene_with_cnmf_on_overlap(tmp_path):
    cellbygene_df, _, _ = _species_cells(overlap=True)

    species_to_gene, cell_by_ecDNA = species_deconvolution.hybrid_deconvolution(
        cellbygene_df, [], n_iter=5, log_dir=str(tmp_path), seed=1
    )

    assigned = set().union(*species_to_gene.values())
    assert assigned == set(cellbygene_df.columns)
    assert _species_of(species_to_gene, "A1") != _species_of(species_to_gene, "B1")
    assert list(cell_by_ecDNA.columns) == [
        f"pred_ecDNA_{i + 1}" for i in range(len(species_to_gene))
    ]
