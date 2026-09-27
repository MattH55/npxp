import numpy as np
import pandas as pd
import pytest

from npi_pharma.interact.score import cosine
from npi_pharma.signatures.build import build_signature, consensus, paired_de


def test_paired_de_recovers_planted_genes():
    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(1000)]
    base = rng.normal(8, 1, (1000, 8))
    post = base + rng.normal(0, 0.2, base.shape)
    post[:10] += 1.5
    post[10:20] -= 1.5
    cols = [f"s{i}" for i in range(8)]
    sig = build_signature(pd.DataFrame(base, genes, cols), pd.DataFrame(post, genes, cols),
                          {"npi_id": "T", "tissue": "adipose", "modality": "diet"})
    assert set(sig.up_genes[:10]) == {f"G{i}" for i in range(10)}
    assert set(sig.down_genes[:10]) == {f"G{i}" for i in range(10, 20)}
    assert sig.quality_flag == "ok" and sig.sample_size == 8


def test_small_n_flagged():
    rng = np.random.default_rng(1)
    genes = [f"G{i}" for i in range(300)]
    a = pd.DataFrame(rng.normal(size=(300, 4)), genes)
    sig = build_signature(a, a + rng.normal(0, 0.1, a.shape), {"npi_id": "T"})
    assert sig.quality_flag == "small_n"


def test_paired_de_needs_two_subjects():
    df = pd.DataFrame({"s": [1.0, 2.0]}, index=["A", "B"])
    with pytest.raises(ValueError):
        paired_de(df, df)


def test_cr_and_exercise_are_two_distinct_ids(fx):
    arms = [s for s in fx["npis"] if s.source_accessions == ["SYNTHETIC_CR_EX_STUDY"]]
    assert sorted(s.sig_id for s in arms) == ["FIXTURE_CR_muscle", "FIXTURE_EX_muscle"]
    assert {s.modality for s in arms} == {"caloric_restriction", "exercise"}
    cr, ex = arms
    assert abs(cosine(cr.z, ex.z)) < 0.3  # "cosine << 1"
    # mouse-case symbols from the builder were normalised onto human symbols
    assert "MTOR" in cr.genes


def test_catalog_keeps_cr_and_exercise_separate():
    from pathlib import Path

    from npi_pharma.ingest.catalog import load_catalog

    entries = load_catalog(Path(__file__).parents[1] / "configs" / "npi_catalog.yaml")
    ids = [e.npi_id for e in entries]
    assert len(ids) == len(set(ids))
    nutr = [e for e in entries if "mdpi.com/2072-6643/15/4/1047" in e.record["source_accessions"][0]]
    assert {e.record["modality"] for e in nutr} == {"caloric_restriction", "exercise"} and len(nutr) == 2
    lcd = next(e for e in entries if e.npi_id == "LCD_adipose_GSE95640")
    assert lcd.quality_flag == "unverified_metadata" and any("not yet curated" in i for i in lcd.issues)


def test_consensus_refuses_mixed_tissues(by_id):
    with pytest.raises(ValueError, match="refusing consensus"):
        consensus([by_id["FIXTURE_LCD_adipose"], by_id["FIXTURE_xenobiotic_diet_liver"]], "bad")
