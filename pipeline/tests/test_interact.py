import numpy as np
import pytest

from npi_pharma.genes import GeneOverlapError, align
from npi_pharma.interact.rank import rank_pairs, select_npis
from npi_pharma.interact.score import score_pair
from npi_pharma.model import PatientState
from npi_pharma.patient.encode import encode_patient


def test_lcd_vs_fixture_drugs_nonzero_and_deterministic(fx, patient, by_id, gene_sets, cfg):
    npi = [by_id["FIXTURE_LCD_adipose"]]
    t1, _ = rank_pairs(patient, npi, fx["drugs"], gene_sets, cfg)
    t2, _ = rank_pairs(patient, npi, list(reversed(fx["drugs"])), gene_sets, cfg)
    assert len(t1) == 50 and t1["composite"].notna().all()
    assert t1["composite"].abs().max() > 0 and t1["composite"].std() > 0
    assert t1["drug_id"].tolist() == t2["drug_id"].tolist()
    assert t1["composite_rank"].tolist() == list(range(1, 51))


def test_fixture_rebuild_is_bitwise_deterministic(fx):
    from npi_pharma.fixtures import build_all

    again = build_all()
    for a, b in zip(fx["drugs"], again["drugs"]):
        assert a.sig_id == b.sig_id and np.array_equal(a.z, b.z)


def test_complementary_patient_gets_high_gain(by_id, gene_sets, cfg):
    n, d = by_id["FIXTURE_LCD_adipose"], by_id["FIXTURE_cmpd_000"]
    genes, v = align({"n": (n.genes, n.z), "d": (d.genes, d.z)})
    s_p = -(v["n"] / np.linalg.norm(v["n"]) + v["d"] / np.linalg.norm(v["d"]))
    s_p += np.random.default_rng(3).normal(0, 0.1 * s_p.std(), len(s_p))
    pat = PatientState("synthetic", "adipose", genes, s_p, s_p, "constructed")
    rep = score_pair(pat, n, d, gene_sets, cfg)
    assert rep["complementarity"] > 0.3
    assert rep["complementarity_gain"] > 0.05
    assert rep["reverse_combo"] > max(rep["reverse_npi"], rep["reverse_drug"])
    # and the pair ranks first against every other fixture drug
    from npi_pharma.fixtures import build_all

    table, _ = rank_pairs(pat, [n], build_all()["drugs"], gene_sets, cfg)
    assert table.iloc[0]["drug_id"] == "FIXTURE_cmpd_000"


def test_pbmc_patient_with_muscle_npi_warns_and_shrinks(fx, by_id, gene_sets, cfg):
    pbmc = encode_patient(fx["patient_expr"], None, "PBMC", reference=fx["reference_expr"])
    rep = score_pair(pbmc, by_id["FIXTURE_EX_muscle"], by_id["FIXTURE_metformin"], gene_sets, cfg)
    assert "tissue_mismatch" in rep["flags"]
    assert any("tissue" in w for w in rep["warnings"])
    assert rep["confidence"] < 1.0
    assert abs(rep["composite"]) < abs(rep["composite_raw"]) or rep["composite_raw"] == 0
    # and rank() drops it by default
    kept, skipped = select_npis(fx["npis"], pbmc, cfg, npi_class="exercise")
    assert kept == [] and skipped


def test_small_n_npi_is_down_weighted(patient, by_id, gene_sets, cfg):
    rep = score_pair(patient, by_id["FIXTURE_psychosocial_blood"], by_id["FIXTURE_metformin"], gene_sets, cfg)
    assert "small_n" in rep["flags"] and "tissue_mismatch" in rep["flags"]
    sh = cfg["shrink"]
    assert "drug_tissue_mismatch" in rep["flags"]  # fixture drugs are cell-line consensus
    assert rep["confidence"] == pytest.approx(sh["small_n"] * sh["tissue_mismatch"] * sh["drug_tissue_mismatch"])


def test_xenobiotic_flag_needs_both_signatures(patient, by_id, gene_sets, cfg):
    npi = by_id["FIXTURE_xenobiotic_diet_liver"]
    both = score_pair(patient, npi, by_id["FIXTURE_pxr_agonist"], gene_sets, cfg)
    only_npi = score_pair(patient, npi, by_id["FIXTURE_rapamycin"], gene_sets, cfg)
    assert "xenobiotic_genes_moved_by_both" in both["flags"]
    assert "xenobiotic_genes_moved_by_both" not in only_npi["flags"]


def test_pathway_antagonism_rapamycin_insulin(patient, by_id, gene_sets, cfg):
    rep = score_pair(patient, by_id["FIXTURE_LCD_adipose"], by_id["FIXTURE_rapamycin"], gene_sets, cfg)
    ins = next(r for r in rep["pathways"] if r["pathway"] == "INSULIN_SIGNALING_CORE")
    assert ins["conflict"] > 0 and ins["joint"] < 0


def test_rank_keeps_unscorable_pairs_visible(patient, by_id, gene_sets, cfg):
    from npi_pharma.model import DRUG, Signature

    tiny = Signature("tiny", DRUG, ["MTOR", "RPTOR"], np.array([1.0, -1.0]), "local_matrix")
    table, _ = rank_pairs(patient, [by_id["FIXTURE_LCD_adipose"]], [tiny, by_id["FIXTURE_metformin"]], gene_sets, cfg)
    row = table[table.drug_id == "tiny"].iloc[0]
    assert np.isnan(row["composite"]) and row["flags"].startswith("error:")
    with pytest.raises(GeneOverlapError):
        score_pair(patient, by_id["FIXTURE_LCD_adipose"], tiny, gene_sets, cfg)


def test_weights_are_exposed(patient, by_id, gene_sets, cfg):
    import copy

    c2 = copy.deepcopy(cfg)
    c2["weights"] = {k: 0.0 for k in cfg["weights"]} | {"pathway_joint": 1.0}
    rep = score_pair(patient, by_id["FIXTURE_LCD_adipose"], by_id["FIXTURE_metformin"], gene_sets, c2)
    assert rep["composite_raw"] == pytest.approx(rep["pathway_joint"])
