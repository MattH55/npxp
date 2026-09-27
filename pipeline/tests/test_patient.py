import numpy as np
import pandas as pd
import pytest

from npi_pharma.patient.encode import encode_patient, reference_offset


def _world(seed=0):
    """Lean and obese groups share a batch; a separate obese cohort holds the patient."""
    rng = np.random.default_rng(seed)
    genes = [f"G{i}" for i in range(400)]
    shift = np.zeros(400)
    shift[:20] = 2.0  # obesity raises G0..G19 by 2 SD
    lean = pd.DataFrame(rng.normal(8, 1, (400, 30)), genes)
    obese = pd.DataFrame(rng.normal(8, 1, (400, 30)) + shift[:, None], genes)
    cohort = pd.DataFrame(rng.normal(8, 1, (400, 40)) + shift[:, None] + 3.0, genes,  # other batch: +3 everywhere
                          columns=[f"p{i}" for i in range(40)])
    return genes, lean, obese, cohort


def test_reference_offset_recovers_effect_size():
    genes, lean, obese, _ = _world()
    d = reference_offset(obese, lean)
    assert d.loc[[f"G{i}" for i in range(20)]].mean() == pytest.approx(2.0, abs=0.3)
    assert abs(d.loc[[f"G{i}" for i in range(20, 400)]].mean()) < 0.1
    with pytest.raises(ValueError, match="samples per group"):
        reference_offset(obese.iloc[:, :3], lean)


def test_offset_places_cohort_patient_against_health_without_batch_shift():
    genes, lean, obese, cohort = _world()
    d = reference_offset(obese, lean)
    p = encode_patient(cohort, "p0", "adipose", offset=d, offset_desc="obese vs lean")
    z = pd.Series(p.disease_vector, index=p.genes)
    assert "health_via_external_offset" in p.flags and "cohort_relative" not in p.flags
    # the +3 batch shift of the cohort does not leak in; the obesity genes stand out
    assert z.loc[[f"G{i}" for i in range(20)]].mean() > 1.5
    assert abs(z.loc[[f"G{i}" for i in range(20, 400)]].mean()) < 0.5
    with pytest.raises(ValueError, match="cohort"):
        encode_patient(cohort, "p0", "adipose", reference=lean, offset=d)
