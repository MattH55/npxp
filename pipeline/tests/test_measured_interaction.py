import numpy as np
import pandas as pd
import pytest

from npi_pharma.interact.measured import (
    KIND_INTERACTION, additivity_summary, build_interaction_signature, interaction_contrast,
)


def _world(seed=0, n_genes=500, n=4, noise=0.2):
    """Four arms where the combination is additive except for planted genes."""
    rng = np.random.default_rng(seed)
    genes = [f"G{i}" for i in range(n_genes)]
    a_eff, b_eff = np.zeros(n_genes), np.zeros(n_genes)
    a_eff[:10] = 1.0                      # agent A moves G0..G9
    b_eff[10:20] = 1.0                    # agent B moves G10..G19
    inter = np.zeros(n_genes)
    inter[:5] = 2.0                       # supra-additive on G0..G4
    inter[20:25] = -2.0                   # sub-additive on G20..G24

    def arm(shift):
        return pd.DataFrame(rng.normal(8, noise, (n_genes, n)) + shift[:, None], index=genes)

    return (arm(np.zeros(n_genes)), arm(a_eff), arm(b_eff), arm(a_eff + b_eff + inter))


def test_contrast_recovers_planted_interaction_and_ignores_additive_genes():
    ctrl, a, b, combo = _world()
    effect, stat = interaction_contrast(ctrl, a, b, combo)
    assert set(effect.nlargest(5).index) == {f"G{i}" for i in range(5)}
    assert set(effect.nsmallest(5).index) == {f"G{i}" for i in range(20, 25)}
    # genes where the combination really is additive sit near zero
    additive = [f"G{i}" for i in range(100, 400)]
    assert effect[additive].abs().mean() < 0.3
    # a purely additive combination leaves no interaction anywhere
    ctrl2, a2, b2, _ = _world(seed=1)
    add_combo = pd.DataFrame(a2.to_numpy() + b2.to_numpy() - ctrl2.to_numpy(),
                             index=a2.index, columns=a2.columns)
    eff2, _ = interaction_contrast(ctrl2, a2, b2, add_combo)
    assert eff2.abs().max() < effect.abs().max() / 2


def test_contrast_is_symmetric_in_the_two_agents():
    ctrl, a, b, combo = _world()
    e1, _ = interaction_contrast(ctrl, a, b, combo)
    e2, _ = interaction_contrast(ctrl, b, a, combo)
    pd.testing.assert_series_equal(e1, e2)


def test_contrast_refuses_too_few_samples_or_genes():
    ctrl, a, b, combo = _world()
    with pytest.raises(ValueError, match=">= 2 samples"):
        interaction_contrast(ctrl.iloc[:, :1], a, b, combo)
    with pytest.raises(ValueError, match="genes shared"):
        interaction_contrast(ctrl.iloc[:50], a.iloc[:50], b.iloc[:50], combo.iloc[:50])


def test_signature_fields_and_small_n_flag():
    ctrl, a, b, combo = _world(n=3)
    meta = {"sig_id": "A + B|GSE1", "agent_a": "A", "agent_b": "B",
            "cell_line": "HCT116", "source_accessions": ["GSE1"]}
    sig = build_interaction_signature(ctrl, a, b, combo, meta)
    assert sig.kind == KIND_INTERACTION and sig.sig_id == "A + B|GSE1"
    assert sig.tissue == "HCT116" and sig.source_accessions == ["GSE1"]
    assert sig.quality_flag == "ok" and sig.sample_size == 3
    assert sig.meta["n_combo"] == 3 and sig.meta["agent_a"] == "A"
    assert set(list(sig.meta["effect"])[:5]) <= {f"G{i}" for i in range(10)}
    # the standardised z keeps the same ordering as the raw effect
    z = sig.as_series()
    assert set(z.nlargest(5).index) == {f"G{i}" for i in range(5)}

    small = build_interaction_signature(ctrl.iloc[:, :2], a.iloc[:, :2], b.iloc[:, :2],
                                        combo.iloc[:, :2], meta)
    assert small.quality_flag == "small_n"


def test_additivity_summary_reports_both_directions():
    ctrl, a, b, combo = _world()
    effect, _ = interaction_contrast(ctrl, a, b, combo)
    s = additivity_summary(effect, top=5)
    assert s["n_genes"] == 500 and s["sd_of_interaction"] > 0
    assert all(v > 0 for v in s["genes_above_additive"].values())
    assert all(v < 0 for v in s["genes_below_additive"].values())
