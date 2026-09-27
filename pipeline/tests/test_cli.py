import json

from npi_pharma.cli import main
from npi_pharma.store import load_patient, load_signatures, save_signatures


def test_store_roundtrip(fx, tmp_path):
    p = save_signatures(fx["npis"], tmp_path / "n.parquet")
    back = {s.sig_id: s for s in load_signatures(p)}
    orig = fx["npis"][0]
    got = back[orig.sig_id]
    assert got.genes == orig.genes and got.tissue == orig.tissue and got.provenance == "synthetic_fixture"
    assert abs(got.z - orig.z).max() < 1e-5


def test_milestone_path_end_to_end(tmp_path):
    out = tmp_path / "demo"
    assert main(["demo", "--out", str(out)]) == 0
    rep = json.loads((out / "report.json").read_text())
    assert rep["evidence_tier"] == "tier_3_mechanism_only"
    s = rep["score"]
    for k in ("complementarity", "orthogonality", "monotherapy_correlation", "pathway_joint", "composite",
              "composite_rank"):
        assert k in s
    assert "SYNTHETIC_FIXTURE" in s["flags"]
    assert rep["explanation"]["pathways_driving"]
    assert (out / "ranked.tsv").read_text().count("\n") == 51
    assert load_patient(out / "patient.npz").reference.startswith("healthy reference")

    assert main(["rank", "--patient", str(out / "patient.npz"), "--npis", str(out / "npis.parquet"),
                 "--drugs", str(out / "drugs.parquet"), "--npi-class", "diet", "--top", "10",
                 "--out", str(out / "rank.tsv"), "--json", str(out / "rank.json")]) == 0
    lines = (out / "rank.tsv").read_text().splitlines()
    assert len(lines) == 11 and all("FIXTURE_LCD_adipose" in ln for ln in lines[1:])  # liver/blood NPIs filtered
    assert json.loads((out / "rank.json").read_text())["top_pairs"]
