import numpy as np
import pytest

from glari import (
    CLINVAR_DB,
    DEFAULT_VARIANTS,
    ModelParams,
    initial_state,
    sample_patient_genomics_mc,
    simulate,
)
from glari.cli import main


def _setup(variants=DEFAULT_VARIANTS, ebv=True, seed=42):
    rng = np.random.default_rng(seed)
    penalties, weights = sample_patient_genomics_mc(list(variants), CLINVAR_DB, rng)
    params = ModelParams(gc_risk_factor=float(weights[0]), ebv_active=ebv)
    return initial_state(penalties[0]), params


def test_sampling_is_reproducible_with_seed():
    a = sample_patient_genomics_mc(list(DEFAULT_VARIANTS), CLINVAR_DB, np.random.default_rng(1), 5)
    b = sample_patient_genomics_mc(list(DEFAULT_VARIANTS), CLINVAR_DB, np.random.default_rng(1), 5)
    np.testing.assert_array_equal(a[0], b[0])


def test_gc_weight_is_product_of_variant_weights():
    _, weights = sample_patient_genomics_mc(list(DEFAULT_VARIANTS), CLINVAR_DB, np.random.default_rng(0))
    assert weights[0] == pytest.approx(2.5 * 1.2)


def test_no_variants_gives_no_penalty():
    penalties, weights = sample_patient_genomics_mc([], CLINVAR_DB, np.random.default_rng(0))
    assert penalties[0] == 0.0 and weights[0] == 1.0


def test_unknown_variant_rejected():
    with pytest.raises(ValueError, match="Unknown ClinVar"):
        sample_patient_genomics_mc(["VCV999999999"], CLINVAR_DB, np.random.default_rng(0))


@pytest.mark.parametrize("variants,ebv", [(DEFAULT_VARIANTS, True), ((), False)])
@pytest.mark.parametrize("treat", [False, True])
def test_simulation_completes_and_stays_in_physical_range(variants, ebv, treat):
    y0, params = _setup(variants, ebv)
    p, f, delta, lam, D = simulate(y0, params, 30.0, treat=treat)
    assert 0.0 <= p <= 1.0
    assert f >= 0.0
    assert params.delta_min <= delta <= params.delta_max
    assert 0.0 <= lam <= 1.0
    assert 0.0 < D <= params.D_max


def test_therapy_lowers_load_and_restores_integrity():
    y0, params = _setup()
    untreated = simulate(y0, params, 30.0, treat=False)
    treated = simulate(y0, params, 30.0, treat=True)
    assert treated[1] < untreated[1]
    assert treated[0] > untreated[0]


def test_therapy_after_horizon_has_no_effect():
    y0, params = _setup()
    np.testing.assert_allclose(simulate(y0, params, 10.0, treat=True), simulate(y0, params, 10.0, treat=False))


def test_cli_runs(capsys):
    main(["--seed", "42"])
    assert "Pipeline Complete" in capsys.readouterr().out


def test_cli_rejects_bad_input():
    with pytest.raises(SystemExit):
        main(["--days", "0"])
    with pytest.raises(SystemExit):
        main(["--variants", "BAD"])
