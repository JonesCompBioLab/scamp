from pathlib import Path

import numpy as np
import pytest
import torch

from scamp.predict import predict_methods

PRETRAINED_MODEL = Path(__file__).resolve().parents[2] / "pretrained_models" / "scamp_model_1.0"
COLO320DM_COPY_NUMBERS = Path(__file__).resolve().parents[2] / "data" / "colo320dm_cna_subset.tsv"


class FakeModel:
    """Stands in for `models.SCAMP`.

    `prepare_copy_numbers` applies the mean filter like the real model and summarizes
    each distribution by its mean and variance; `proba` turns that summary into an
    ecDNA probability with `probability_fn`. Every call to `prepare_copy_numbers` is
    recorded in `calls`.

    The probability function is a dummy function that merely returns a fixed
    probability.
    """

    def __init__(self, probability_fn=None):
        self.probability_fn = probability_fn or (lambda mean, var: 0.9)
        self.calls = []

    def prepare_copy_numbers(
        self, copy_numbers, genes, min_copy_number, max_percentile, filter_copy_number=3
    ):
        self.calls.append(
            {
                "copy_numbers": np.array(copy_numbers, dtype=float),
                "min_copy_number": min_copy_number,
                "max_percentile": max_percentile,
                "filter_copy_number": filter_copy_number,
            }
        )
        means = np.mean(copy_numbers, axis=0)
        kii = np.where(means >= filter_copy_number)[0]
        X = np.array(
            [[np.mean(copy_numbers[:, i]), np.var(copy_numbers[:, i])] for i in kii]
        ).reshape(len(kii), 2)
        return X, np.asarray(genes)[kii]

    def proba(self, X):
        p = self.probability_fn(float(X[0, 0]), float(X[0, 1]))
        return torch.tensor([[1 - p, p]])


def _high_variance_is_ecdna(mean, var):
    return 0.9 if var > 5 else 0.1


def _diploid_cells(n=1000, seed=0):
    """Cells with a narrow copy-number distribution around 3."""
    return np.random.default_rng(seed).normal(3, 0.2, size=n)


def _diploid_and_amplified_cells(seed=0):
    """800 cells near copy-number 3, and 200 with a wide distribution around 20."""
    rng = np.random.default_rng(seed)
    return np.concatenate([rng.normal(3, 0.2, size=800), rng.normal(20, 4, size=200)])


GENES = np.array(["MYC"])


# _vector_seed


def test_vector_seed_is_deterministic_and_ignores_shape():
    values = [2.0, 3.5, 10.25]

    assert predict_methods._vector_seed(values) == predict_methods._vector_seed(values)
    assert predict_methods._vector_seed(values) == predict_methods._vector_seed(
        np.array(values).reshape(-1, 1)
    )


def test_vector_seed_rounds_to_three_decimals_and_treats_negative_zero_as_zero():
    assert predict_methods._vector_seed([1.0001, 2.0]) == predict_methods._vector_seed(
        [1.0, 2.0]
    )
    assert predict_methods._vector_seed([-0.0, 1.0]) == predict_methods._vector_seed(
        [0.0, 1.0]
    )


def test_vector_seed_differs_for_different_vectors():
    seed = predict_methods._vector_seed([1.0, 2.0])

    assert seed != predict_methods._vector_seed([1.0, 2.01])
    assert 0 <= seed < 2**64


# _trim_cells
def test_trim_cells_keeps_all_cells_by_default():
    result = predict_methods._trim_cells(np.array([[3.0], [1.0], [2.0]]))

    assert result.ndim == 1
    np.testing.assert_array_equal(result, [3.0, 1.0, 2.0])


def test_trim_cells_drops_cells_above_max_percentile():
    # The 90th percentile of 1..10 is 9.1
    result = predict_methods._trim_cells(np.arange(1, 11), max_percentile=90)

    np.testing.assert_array_equal(result, np.arange(1, 10))


def test_trim_cells_min_copy_number_is_strict():
    result = predict_methods._trim_cells([1.0, 2.0, 2.5, 4.0], min_copy_number=2.0)

    np.testing.assert_array_equal(result, [2.5, 4.0])


def test_trim_cells_takes_the_percentile_over_all_cells_before_the_minimum():
    # The 50th percentile is taken over all six cells (3.5), not only those above 2
    result = predict_methods._trim_cells(
        [1, 2, 3, 4, 5, 6], max_percentile=50, min_copy_number=2
    )

    np.testing.assert_array_equal(result, [3])


def test_trim_cells_returns_none_when_kept_mean_is_below_filter():
    assert predict_methods._trim_cells([2.1, 2.2, 2.3], filter_copy_number=2.5) is None


def test_trim_cells_keeps_cells_whose_mean_equals_the_filter():
    result = predict_methods._trim_cells([2.0, 3.0], filter_copy_number=2.5)

    np.testing.assert_array_equal(result, [2.0, 3.0])


def test_trim_cells_returns_empty_array_not_none_when_no_cells_are_left():
    result = predict_methods._trim_cells(
        [1.0, 1.5], min_copy_number=2.0, filter_copy_number=2.5
    )

    assert result is not None
    assert len(result) == 0


# Test GMM fitting
def test_gmm_fit_uses_given_number_of_components():
    gmm = predict_methods._GMM_fit(_diploid_and_amplified_cells().reshape(-1, 1), n_components=3)

    assert gmm.n_components == 3


def test_gmm_fit_chooses_one_component_for_unimodal_data():
    gmm = predict_methods._GMM_fit(
        np.random.default_rng(0).normal(10, 2, size=1000).reshape(-1, 1)
    )

    assert gmm.n_components == 1


def test_gmm_fit_chooses_one_component_when_bic_is_negative():
    # A very narrow distribution has a large log-likelihood, so its BIC is negative
    cells = np.random.default_rng(0).normal(3, 0.01, size=1000).reshape(-1, 1)
    bic = predict_methods.GaussianMixture(n_components=1, random_state=0).fit(cells).bic(cells)
    assert bic < 0

    gmm = predict_methods._GMM_fit(cells)

    assert gmm.n_components == 1


def test_gmm_fit_chooses_several_components_when_the_best_bic_is_negative():
    # Two narrow modes: one distribution fits badly (positive BIC), two fit very well
    # (negative BIC). The old bic[0] / min(bic) ratio was negative here and chose 1.
    rng = np.random.default_rng(0)
    cells = np.concatenate([rng.normal(3, 0.01, 500), rng.normal(6, 0.01, 500)]).reshape(-1, 1)
    bics = [
        predict_methods.GaussianMixture(n_components=k, random_state=0).fit(cells).bic(cells)
        for k in (1, 2)
    ]
    assert bics[0] > 0 > bics[1]

    gmm = predict_methods._GMM_fit(cells)

    assert gmm.n_components > 1


def test_gmm_fit_chooses_several_components_for_separated_modes():
    gmm = predict_methods._GMM_fit(_diploid_and_amplified_cells().reshape(-1, 1))

    assert gmm.n_components > 1


# Predict with SCAMP model
def test_scamp_predict_scores_a_distribution_that_passes_the_filter():
    model = FakeModel(lambda mean, var: 0.8)

    prediction, probability = predict_methods.SCAMP_ecDNA_predict(
        model, [3.0, 4.0, 5.0], [], GENES, decision_rule=0.5
    )

    assert prediction
    assert probability == pytest.approx(0.8)
    assert model.calls[0]["copy_numbers"].shape == (3, 1)


def test_scamp_predict_decision_rule_is_strict():
    model = FakeModel(lambda mean, var: 0.5)

    prediction, _ = predict_methods.SCAMP_ecDNA_predict(
        model, [3.0, 4.0], [], GENES, decision_rule=0.5
    )

    assert not prediction


def test_scamp_predict_passes_preprocessing_to_the_model():
    model = FakeModel()

    predict_methods.SCAMP_ecDNA_predict(
        model,
        [3.0, 4.0],
        [],
        GENES,
        min_copy_number=1.0,
        max_percentile=95.0,
        filter_copy_number=0.5,
    )

    call = model.calls[0]
    assert (call["min_copy_number"], call["max_percentile"], call["filter_copy_number"]) == (
        1.0,
        95.0,
        0.5,
    )


def test_scamp_predict_returns_false_and_logs_when_filtered_out():
    out_log = []

    result = predict_methods.SCAMP_ecDNA_predict(
        FakeModel(), [2.0, 2.1], out_log, GENES, verbose=True, filter_copy_number=2.5
    )

    assert result == (False, 0)
    assert out_log == [" MYC did not pass filter"]


def test_scamp_predict_does_not_log_when_not_verbose():
    out_log = []

    predict_methods.SCAMP_ecDNA_predict(
        FakeModel(), [2.0, 2.1], out_log, GENES, filter_copy_number=2.5
    )

    assert out_log == []


# Test GMM ecDNA prediction
def test_gmm_predict_calls_the_wide_amplified_distribution_ecdna():
    model = FakeModel(_high_variance_is_ecdna)

    prediction, probability, mean, var, weight = predict_methods.GMM_ecDNA_predict(
        model, _diploid_and_amplified_cells(), [], genes=GENES, trim_before_gmm=False
    )

    assert prediction
    assert probability == pytest.approx(0.9)
    assert mean > 10
    assert var > 5
    assert 0 < weight < 1


def test_gmm_predict_negative_reports_highest_probability_checked():
    model = FakeModel(lambda mean, var: 0.2)

    result = predict_methods.GMM_ecDNA_predict(
        model, _diploid_and_amplified_cells(), [], genes=GENES, trim_before_gmm=False
    )

    assert result[0] is False
    assert result[1] == pytest.approx(0.2)
    assert result[2:] == (-1, -1, -1)


def test_gmm_predict_decision_rule_is_inclusive():
    model = FakeModel(lambda mean, var: 0.5)

    prediction, *_ = predict_methods.GMM_ecDNA_predict(
        model,
        _diploid_cells(),
        [],
        genes=GENES,
        decision_rule=0.5,
        trim_before_gmm=False,
    )

    assert prediction


def test_gmm_predict_scores_sampled_distributions_without_model_preprocessing():
    model = FakeModel(lambda mean, var: 0.1)

    predict_methods.GMM_ecDNA_predict(
        model, _diploid_and_amplified_cells(), [], genes=GENES, trim_before_gmm=False
    )

    assert model.calls
    for call in model.calls:
        assert call["copy_numbers"].shape == (1000, 1)
        assert call["min_copy_number"] == -np.inf
        assert call["max_percentile"] == 100
        assert call["filter_copy_number"] == -np.inf


def test_gmm_predict_is_reproducible():
    cells = _diploid_and_amplified_cells()

    first = predict_methods.GMM_ecDNA_predict(
        FakeModel(_high_variance_is_ecdna), cells, [], genes=GENES, trim_before_gmm=False
    )
    second = predict_methods.GMM_ecDNA_predict(
        FakeModel(_high_variance_is_ecdna), cells, [], genes=GENES, trim_before_gmm=False
    )

    assert first == second


def test_gmm_predict_scores_each_distribution_n_draws_times():
    model = FakeModel(lambda mean, var: 0.1)

    predict_methods.GMM_ecDNA_predict(
        model,
        _diploid_cells(),
        [],
        genes=GENES,
        n_components=1,
        n_draws=3,
        trim_before_gmm=False,
    )

    assert len(model.calls) == 3


def test_gmm_predict_skips_distributions_below_min_weight():
    model = FakeModel(lambda mean, var: 0.1)

    # 800 / 200 cells, so only the larger distribution has weight above 0.5
    predict_methods.GMM_ecDNA_predict(
        model,
        _diploid_and_amplified_cells(),
        [],
        genes=GENES,
        n_components=2,
        min_weight=0.5,
        trim_before_gmm=False,
    )

    assert len(model.calls) == 1


def test_gmm_predict_not_called_when_trimmed_mean_is_below_filter():
    out_log = []

    result = predict_methods.GMM_ecDNA_predict(
        FakeModel(),
        np.full(100, 2.2),
        out_log,
        genes=GENES,
        min_copy_number=2.0,
        filter_copy_number=2.5,
    )

    assert result == (False, 0.0, -1, -1, -1)
    assert out_log == ["MYC: mean copy-number after GMM trimming is below 2.5, not called"]


def test_gmm_predict_not_called_when_too_few_cells_are_left_after_trimming():
    out_log = []
    cells = np.concatenate([np.full(100, 1.0), np.full(5, 10.0)])

    result = predict_methods.GMM_ecDNA_predict(
        FakeModel(), cells, out_log, genes=GENES, min_copy_number=2.0, max_percentile=100
    )

    assert result == (False, 0.0, -1, -1, -1)
    assert out_log == ["MYC: only 5 cells left after GMM trimming, not called"]


def test_gmm_predict_single_distribution_uses_scamp_scoring_on_all_cells():
    model = FakeModel(lambda mean, var: 0.7)
    cells = _diploid_cells()

    prediction, probability, mean, var, weight = predict_methods.GMM_ecDNA_predict(
        model, cells, [], genes=GENES, single_dist_scamp=True, trim_before_gmm=False
    )

    assert len(model.calls) == 1
    call = model.calls[0]
    np.testing.assert_array_equal(call["copy_numbers"].ravel(), cells)
    assert (call["min_copy_number"], call["max_percentile"], call["filter_copy_number"]) == (
        2.0,
        99.9,
        2.5,
    )
    assert prediction
    assert probability == pytest.approx(0.7)
    assert mean == pytest.approx(3, abs=0.1)
    assert weight == pytest.approx(1.0)


# KNN_ecDNA_predict
def test_knn_predict_counts_every_cell_when_all_neighbourhoods_are_ecdna():
    cells = _diploid_and_amplified_cells()

    prediction, num_positive = predict_methods.KNN_ecDNA_predict(
        FakeModel(lambda mean, var: 0.9),
        cells,
        [],
        GENES,
        trim_before_gmm=False,
        knn_min_mean=-np.inf,
        knn_min_var=-np.inf,
    )

    assert prediction
    assert num_positive == len(cells)


def test_knn_predict_gate_skips_neighbourhoods():
    model = FakeModel()

    result = predict_methods.KNN_ecDNA_predict(
        model,
        _diploid_and_amplified_cells(),
        [],
        GENES,
        trim_before_gmm=False,
        knn_min_mean=1000,
    )

    assert result == (False, 0)
    assert model.calls == []


def test_knn_predict_neighbourhoods_have_at_most_k_plus_one_cells():
    model = FakeModel(lambda mean, var: 0.1)

    predict_methods.KNN_ecDNA_predict(
        model,
        _diploid_cells(500),
        [],
        GENES,
        n_components=1,
        k=20,
        trim_before_gmm=False,
        knn_min_mean=-np.inf,
        knn_min_var=-np.inf,
    )

    sizes = {len(call["copy_numbers"]) for call in model.calls}
    assert max(sizes) == 21


def test_knn_predict_scores_neighbourhoods_without_model_preprocessing():
    model = FakeModel(lambda mean, var: 0.1)

    predict_methods.KNN_ecDNA_predict(
        model,
        _diploid_cells(200),
        [],
        GENES,
        n_components=1,
        trim_before_gmm=False,
        knn_min_mean=-np.inf,
        knn_min_var=-np.inf,
    )

    assert model.calls
    for call in model.calls:
        assert call["min_copy_number"] == -np.inf
        assert call["max_percentile"] == 100
        assert call["filter_copy_number"] == -np.inf


def test_knn_predict_subsamples_to_4000_cells():
    cells = np.random.default_rng(0).normal(10, 2, size=6000)

    _, num_positive = predict_methods.KNN_ecDNA_predict(
        FakeModel(lambda mean, var: 0.9),
        cells,
        [],
        GENES,
        n_components=1,
        trim_before_gmm=False,
        knn_min_mean=-np.inf,
        knn_min_var=-np.inf,
    )

    assert num_positive == 4000


def test_knn_predict_needs_more_than_ecdna_percentage_thresh_positive():
    # Only neighbourhoods with high variance count as ecDNA; none do here
    prediction, num_positive = predict_methods.KNN_ecDNA_predict(
        FakeModel(_high_variance_is_ecdna),
        _diploid_cells(),
        [],
        GENES,
        n_components=1,
        trim_before_gmm=False,
        knn_min_mean=-np.inf,
        knn_min_var=-np.inf,
        ecDNA_percentage_thresh=0.001,
    )

    assert not prediction
    assert num_positive == 0


def test_knn_predict_is_reproducible():
    cells = np.random.default_rng(1).normal(10, 3, size=5000)
    kwargs = dict(n_components=2, trim_before_gmm=False, knn_min_mean=-np.inf, knn_min_var=-np.inf)

    first = predict_methods.KNN_ecDNA_predict(
        FakeModel(_high_variance_is_ecdna), cells, [], GENES, **kwargs
    )
    second = predict_methods.KNN_ecDNA_predict(
        FakeModel(_high_variance_is_ecdna), cells, [], GENES, **kwargs
    )

    assert first == second


def test_knn_predict_not_called_when_trimmed_mean_is_below_filter():
    out_log = []

    result = predict_methods.KNN_ecDNA_predict(
        FakeModel(), np.full(100, 2.2), out_log, GENES, filter_copy_number=2.5
    )

    assert result == (False, 0)
    assert out_log == ["MYC: mean copy-number after GMM trimming is below 2.5, not called"]


def test_knn_predict_not_called_when_too_few_cells_are_left_after_trimming():
    out_log = []
    cells = np.concatenate([np.full(100, 1.0), np.full(5, 10.0)])

    result = predict_methods.KNN_ecDNA_predict(
        FakeModel(), cells, out_log, GENES, max_percentile=100
    )

    assert result == (False, 0)
    assert out_log == ["MYC: only 5 cells left after GMM trimming, not called"]


def test_knn_predict_single_distribution_uses_scamp_scoring_on_all_cells():
    model = FakeModel(lambda mean, var: 0.7)
    cells = _diploid_cells()

    result = predict_methods.KNN_ecDNA_predict(
        model, cells, [], GENES, single_dist_scamp=True, trim_before_gmm=False
    )

    assert result == (True, -1)
    assert len(model.calls) == 1
    call = model.calls[0]
    np.testing.assert_array_equal(call["copy_numbers"].ravel(), cells)
    assert (call["min_copy_number"], call["max_percentile"], call["filter_copy_number"]) == (
        2.0,
        99.9,
        2.5,
    )

# Test Pretrained model

@pytest.mark.optional
@pytest.mark.skipif(not PRETRAINED_MODEL.exists(), reason="pretrained model not available")
def test_scamp_predict_with_the_pretrained_model_returns_a_probability():
    from scamp import models

    model = models.SCAMP.load(str(PRETRAINED_MODEL))

    prediction, probability = predict_methods.SCAMP_ecDNA_predict(
        model, _diploid_and_amplified_cells(), [], GENES
    )

    assert isinstance(prediction, bool)
    assert 0.0 <= probability <= 1.0


@pytest.fixture(scope="module")
def colo320dm_and_pretrained_model():
    from scamp import io, models

    if not (PRETRAINED_MODEL.exists() and COLO320DM_COPY_NUMBERS.exists()):
        pytest.skip("pretrained model or COLO320DM copy-numbers not available")
    counts_df = io.read_copy_numbers_file(str(COLO320DM_COPY_NUMBERS))
    model = models.SCAMP.load(str(PRETRAINED_MODEL))
    return counts_df, model


@pytest.mark.optional
@pytest.mark.parametrize("gene", ["MYC", "PVT1"])
def test_pretrained_model_calls_colo320dm_myc_and_pvt1_ecdna(
    colo320dm_and_pretrained_model, gene
):
    # MYC and PVT1 are on ecDNA in COLO320DM
    counts_df, model = colo320dm_and_pretrained_model

    prediction, probability = predict_methods.SCAMP_ecDNA_predict(
        model, counts_df[gene], [], np.array([gene]), decision_rule=0.5
    )

    assert probability > 0.5
    assert prediction


@pytest.mark.optional
@pytest.mark.parametrize("gene", ["KLHL29", "ASAP1"])
def test_pretrained_model_does_not_call_colo320dm_controls_ecdna(
    colo320dm_and_pretrained_model, gene
):
    # Elevated copy-number in this subset (mean ~11 and ~8) but outside the MYC/PVT1
    # ecDNA, so they should not be called
    counts_df, model = colo320dm_and_pretrained_model

    prediction, probability = predict_methods.SCAMP_ecDNA_predict(
        model, counts_df[gene], [], np.array([gene]), decision_rule=0.5
    )

    assert probability < 0.5
    assert not prediction
