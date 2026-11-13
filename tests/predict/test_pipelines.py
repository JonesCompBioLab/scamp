from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import torch

from scamp.predict import pipelines

PRETRAINED_MODEL = Path(__file__).resolve().parents[2] / "pretrained_models" / "scamp_model_1.0"
COLO320DM_COPY_NUMBERS = Path(__file__).resolve().parents[2] / "data" / "colo320dm_cna_subset.tsv"

# Genes on the MYC/PVT1 ecDNA in COLO320DM, and elevated genes off of it
COLO320DM_ECDNA_GENES = ["MYC", "PVT1"]
COLO320DM_CONTROL_GENES = ["KLHL29", "ASAP1"]


class FakeModel:
    """Stands in for `models.SCAMP` so `predict` can run without the pretrained model.

    `prepare_copy_numbers` applies the mean filter like the real model and summarizes
    each distribution by its mean and variance; `proba` calls a gene ecDNA when its
    variance is above 5. Every call to `prepare_copy_numbers` is recorded in `calls`.
    """

    def __init__(self):
        self.calls = []

    def prepare_copy_numbers(
        self, copy_numbers, genes, min_copy_number, max_percentile, filter_copy_number=3
    ):
        self.calls.append({"genes": list(genes), "copy_numbers": np.array(copy_numbers)})
        means = np.mean(copy_numbers, axis=0)
        kii = np.where(means >= filter_copy_number)[0]
        X = np.array(
            [[np.mean(copy_numbers[:, i]), np.var(copy_numbers[:, i])] for i in kii]
        ).reshape(len(kii), 2)
        return X, np.asarray(genes)[kii]

    def proba(self, X):
        p = 0.9 if float(X[0, 1]) > 5 else 0.1
        return torch.tensor([[1 - p, p]])


@pytest.fixture
def fake_model(monkeypatch):
    """Makes `predict` load a `FakeModel` instead of reading a saved model."""
    model = FakeModel()
    monkeypatch.setattr(pipelines.models.SCAMP, "load", lambda path: model)
    return model


def _counts(n_cells=500, seed=0):
    """Cell-by-gene copy-numbers with known calls under `FakeModel`.

    AMP is wide around 20 (ecDNA), AMP_COPY is identical to AMP, FLAT is narrow
    around 3 (not ecDNA) and LOW has mean 1 (below the default filter).
    """
    rng = np.random.default_rng(seed)
    amp = rng.normal(20, 5, n_cells)
    return pd.DataFrame(
        {
            "AMP": amp,
            "AMP_COPY": amp,
            "FLAT": rng.normal(3, 0.2, n_cells),
            "LOW": rng.normal(1, 0.2, n_cells),
        },
        index=[f"cell{i}" for i in range(n_cells)],
    )


def _predict(counts_df, predict_method="SCAMP", **kwargs):
    kwargs = {"decision_rule": 0.5, "filter_copy_number": 2.5, **kwargs}
    return pipelines.predict(
        counts_df, [], "unused", None, predict_method=predict_method, **kwargs
    )


def _cluster(predict_df, counts_df, deconvolution_method="hier", **kwargs):
    """Runs `cluster` with the CLI's defaults, overridden by `kwargs`."""
    args = {
        "filename": "test",
        "hier_ddist": 1.3,
        "overlap_thresh": 0.55,
        "cNMF_num_iter": 5,
        "error_w": 0.25,
        "score_cutoff": 3,
        "ecDNA_remove_cutoff": 0.01,
        "density_threshold": 0.1,
        "max_species": None,
        "num_ecDNA": None,
        "cNMF_log_dir": None,
        "out_log": [],
        "verbose": False,
        "seed": 1,
        **kwargs,
    }
    return pipelines.cluster(
        predict_df=predict_df,
        counts_df=counts_df,
        deconvolution_method=deconvolution_method,
        **args,
    )


def _two_species_counts(n_cells=300, seed=0):
    """Copy-numbers of two independent ecDNA species, A (A1-A3) and B (B1-B3)."""
    rng = np.random.default_rng(seed)
    a = rng.gamma(2.0, 5.0, n_cells)
    b = rng.gamma(2.0, 5.0, n_cells)
    data = {}
    for i, copies in enumerate([1, 1, 2]):
        data[f"A{i + 1}"] = 2 + copies * a + rng.normal(0, 0.05, n_cells)
        data[f"B{i + 1}"] = 2 + copies * b + rng.normal(0, 0.05, n_cells)
    return pd.DataFrame(data, index=[f"cell{i}" for i in range(n_cells)]).clip(lower=2)


def _all_positive(counts_df):
    return pd.DataFrame({"gene": counts_df.columns, "prediction": True})


def _species_of(predict_df):
    return dict(zip(predict_df["gene"], predict_df["Species"]))


# predict


def test_predict_drops_genes_at_or_below_the_filter(fake_model):
    prediction_df, _ = _predict(_counts())

    assert prediction_df["gene"].tolist() == ["AMP", "AMP_COPY", "FLAT"]


def test_predict_calls_genes_with_the_model(fake_model):
    prediction_df, _ = _predict(_counts())

    calls = dict(zip(prediction_df["gene"], prediction_df["prediction"]))
    assert calls == {"AMP": True, "AMP_COPY": True, "FLAT": False}
    assert prediction_df.set_index("gene").loc["AMP", "Probability"] == pytest.approx(0.9)


def test_predict_reuses_the_result_of_identical_copy_numbers(fake_model):
    out_log = []
    pipelines.predict(_counts(), out_log, "unused", None, 0.5, 2.5, "SCAMP")

    # AMP_COPY is never scored, and the log names the gene whose result it reuses
    assert [call["genes"] for call in fake_model.calls] == [["AMP"], ["FLAT"]]
    assert any(line.startswith("Skipping AMP_COPY") and "'AMP'" in line for line in out_log)


def test_predict_reports_the_mean_and_variance_of_each_gene(fake_model):
    counts_df = _counts()
    prediction_df, _ = _predict(counts_df)

    expected = counts_df[["AMP", "AMP_COPY", "FLAT"]]
    np.testing.assert_allclose(prediction_df["mean"], expected.mean())
    np.testing.assert_allclose(prediction_df["var"], expected.var(ddof=0))


@pytest.mark.parametrize(
    "predict_method, columns",
    [
        ("SCAMP", ["gene", "mean", "var", "Probability", "prediction"]),
        ("GMM", ["gene", "mean", "var", "Probability", "Weight", "prediction"]),
        ("KNN", ["gene", "mean", "var", "NumberOfPositiveCells", "prediction"]),
    ],
)
def test_predict_columns_depend_on_the_method(fake_model, predict_method, columns):
    prediction_df, _ = _predict(_counts()[["AMP", "FLAT"]], predict_method=predict_method)

    assert prediction_df.columns.tolist() == columns


def test_predict_keeps_only_whitelisted_cells(fake_model, tmp_path):
    counts_df = _counts()
    whitelist = counts_df.index[:100]
    whitelist_file = tmp_path / "whitelist.txt"
    pd.Series(list(whitelist) + ["not_a_cell"]).to_csv(whitelist_file, index=False, header=False)

    prediction_df, _ = pipelines.predict(
        counts_df, [], "unused", str(whitelist_file), 0.5, 2.5, "SCAMP"
    )

    # The model and the reported statistics only see the whitelisted cells
    assert all(len(call["copy_numbers"]) == 100 for call in fake_model.calls)
    np.testing.assert_allclose(
        prediction_df["mean"], counts_df.loc[whitelist, ["AMP", "AMP_COPY", "FLAT"]].mean()
    )


def test_predict_raises_on_an_unknown_method(fake_model):
    with pytest.raises(ValueError, match="Unknown predict method"):
        _predict(_counts(), predict_method="NOT_A_METHOD")


def test_predict_raises_when_no_gene_passes_the_filter(fake_model):
    with pytest.raises(ValueError, match="No genes had a mean copy-number above 100"):
        _predict(_counts(), filter_copy_number=100)


# cluster


def test_cluster_without_ecdna_genes_assigns_no_species():
    counts_df = _counts()
    predict_df = pd.DataFrame({"gene": ["AMP", "FLAT"], "prediction": [False, False]})

    predict_df, cell_by_ecDNA, _ = _cluster(predict_df, counts_df)

    assert predict_df["Species"].tolist() == ["None", "None"]
    assert cell_by_ecDNA.shape == (len(counts_df), 0)
    assert cell_by_ecDNA.index.equals(counts_df.index)


def test_cluster_single_profile_is_one_species():
    counts_df = _counts()
    predict_df = pd.DataFrame(
        {"gene": ["AMP", "AMP_COPY", "FLAT"], "prediction": [True, True, False]}
    )

    predict_df, cell_by_ecDNA, _ = _cluster(predict_df, counts_df)

    # AMP_COPY collapses onto AMP, and the species amount is the copy-number above 2
    assert _species_of(predict_df) == {
        "AMP": "pred_ecDNA_1",
        "AMP_COPY": "pred_ecDNA_1",
        "FLAT": "None",
    }
    assert cell_by_ecDNA.columns.tolist() == ["pred_ecDNA_1"]
    np.testing.assert_allclose(
        cell_by_ecDNA["pred_ecDNA_1"], (counts_df["AMP"] - 2).clip(lower=0)
    )


@pytest.mark.parametrize("deconvolution_method", ["hier", "auto"])
def test_cluster_separates_two_species(deconvolution_method):
    counts_df = _two_species_counts()

    predict_df, cell_by_ecDNA, _ = _cluster(
        _all_positive(counts_df), counts_df, deconvolution_method
    )

    species = _species_of(predict_df)
    assert len({species["A1"], species["A2"], species["A3"]}) == 1
    assert len({species["B1"], species["B2"], species["B3"]}) == 1
    assert species["A1"] != species["B1"]
    assert sorted(cell_by_ecDNA.columns) == sorted({species["A1"], species["B1"]})


def test_cluster_gives_duplicate_genes_the_species_of_the_kept_gene():
    counts_df = _two_species_counts()
    counts_df["A1_COPY"] = counts_df["A1"]

    predict_df, _, _ = _cluster(_all_positive(counts_df), counts_df)

    species = _species_of(predict_df)
    assert species["A1_COPY"] == species["A1"]


def test_cluster_with_cnmf_separates_two_species(tmp_path):
    counts_df = _two_species_counts()

    predict_df, cell_by_ecDNA, _ = _cluster(
        _all_positive(counts_df), counts_df, "cNMF", num_ecDNA=2, cNMF_log_dir=str(tmp_path)
    )

    species = _species_of(predict_df)
    assert species["A1"] == species["A2"] == species["A3"]
    assert species["B1"] == species["B2"] == species["B3"]
    assert species["A1"] != species["B1"]
    assert sorted(cell_by_ecDNA.columns) == ["pred_ecDNA_1", "pred_ecDNA_2"]


def test_cluster_raises_on_an_unknown_method():
    counts_df = _two_species_counts()

    with pytest.raises(ValueError, match="Unknown deconvolution method"):
        _cluster(_all_positive(counts_df), counts_df, "NOT_A_METHOD")


# Pretrained model on COLO320DM
@pytest.fixture(scope="module")
def colo320dm_copy_numbers():
    from scamp import io

    if not (PRETRAINED_MODEL.exists() and COLO320DM_COPY_NUMBERS.exists()):
        pytest.skip("pretrained model or COLO320DM copy-numbers not available")
    return io.read_copy_numbers_file(str(COLO320DM_COPY_NUMBERS))


@pytest.fixture(scope="module")
def colo320dm_scamp_predictions(colo320dm_copy_numbers):
    prediction_df, _ = pipelines.predict(
        colo320dm_copy_numbers, [], str(PRETRAINED_MODEL), None, 0.5, 2.5, "SCAMP"
    )
    return prediction_df


@pytest.mark.optional
def test_pretrained_scamp_predicts_every_colo320dm_gene(
    colo320dm_copy_numbers, colo320dm_scamp_predictions
):
    # Every gene in the subset has a mean copy-number above 2.5
    assert colo320dm_scamp_predictions["gene"].tolist() == colo320dm_copy_numbers.columns.tolist()
    assert colo320dm_scamp_predictions["Probability"].between(0, 1).all()


@pytest.mark.optional
def test_pretrained_scamp_calls_colo320dm_ecdna_and_not_controls(colo320dm_scamp_predictions):
    predictions = colo320dm_scamp_predictions.set_index("gene")

    for gene in COLO320DM_ECDNA_GENES:
        assert predictions.loc[gene, "prediction"]
        assert predictions.loc[gene, "Probability"] > 0.5
    for gene in COLO320DM_CONTROL_GENES:
        assert not predictions.loc[gene, "prediction"]
        assert predictions.loc[gene, "Probability"] < 0.5


@pytest.mark.optional
@pytest.mark.parametrize("predict_method", ["GMM", "KNN"])
def test_pretrained_gmm_and_knn_call_colo320dm_ecdna_and_not_controls(
    colo320dm_copy_numbers, predict_method
):
    # Uses the CLI's defaults for scoring and trimming
    counts_df = colo320dm_copy_numbers[COLO320DM_ECDNA_GENES + COLO320DM_CONTROL_GENES]

    prediction_df, _ = pipelines.predict(
        counts_df,
        [],
        str(PRETRAINED_MODEL),
        None,
        0.5,
        2.5,
        predict_method,
        trim_before_gmm=True,
        single_dist_scamp=True,
    )

    calls = dict(zip(prediction_df["gene"], prediction_df["prediction"]))
    assert calls == {"MYC": True, "PVT1": True, "KLHL29": False, "ASAP1": False}


@pytest.mark.optional
def test_pretrained_scamp_with_a_whitelist_uses_only_listed_cells(
    colo320dm_copy_numbers, tmp_path
):
    whitelist = colo320dm_copy_numbers.index[:5000]
    whitelist_file = tmp_path / "whitelist.txt"
    pd.Series(whitelist).to_csv(whitelist_file, index=False, header=False)
    counts_df = colo320dm_copy_numbers[COLO320DM_ECDNA_GENES + COLO320DM_CONTROL_GENES]

    prediction_df, _ = pipelines.predict(
        counts_df, [], str(PRETRAINED_MODEL), str(whitelist_file), 0.5, 2.5, "SCAMP"
    )

    np.testing.assert_allclose(prediction_df["mean"], counts_df.loc[whitelist].mean())
    calls = dict(zip(prediction_df["gene"], prediction_df["prediction"]))
    assert calls == {"MYC": True, "PVT1": True, "KLHL29": False, "ASAP1": False}


@pytest.mark.optional
def test_pretrained_scamp_and_hierarchical_put_the_colo320dm_amplicon_on_one_species(
    colo320dm_copy_numbers, colo320dm_scamp_predictions
):
    # MYC, PVT1 and the genes between them are amplified on a single ecDNA in COLO320DM
    predict_df, cell_by_ecDNA, _ = _cluster(
        colo320dm_scamp_predictions.copy(), colo320dm_copy_numbers, "hier"
    )

    species = _species_of(predict_df)
    ecdna_genes = predict_df.loc[predict_df["prediction"], "gene"]
    assert {species[gene] for gene in ecdna_genes} == {"pred_ecDNA_1"}
    assert all(species[gene] == "None" for gene in COLO320DM_CONTROL_GENES)
    assert cell_by_ecDNA.columns.tolist() == ["pred_ecDNA_1"]
    assert cell_by_ecDNA.index.equals(colo320dm_copy_numbers.index)
    assert (cell_by_ecDNA["pred_ecDNA_1"] >= 0).all()
