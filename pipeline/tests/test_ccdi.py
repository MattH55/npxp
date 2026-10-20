import hashlib
import io
import json

import numpy as np
import pandas as pd
import pytest

from npi_pharma.ccdi import annotate, cohorts, manifest, openpedcan


# ---------------------------------------------------------------- fetch / manifest

class _Resp(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def _fake_bucket(monkeypatch, files: dict[str, bytes], md5_override: dict[str, str] | None = None):
    md5 = {k: hashlib.md5(v).hexdigest() for k, v in files.items()} | (md5_override or {})
    pages = {f"{openpedcan.BUCKET}/vX/md5sum.txt": "".join(f"{h}  {k}\n" for k, h in md5.items()).encode()}
    pages |= {f"{openpedcan.BUCKET}/vX/{k}": v for k, v in files.items()}

    def urlopen(url, timeout=None):
        if url not in pages:
            raise OSError(f"404 {url}")
        return _Resp(pages[url])

    monkeypatch.setattr(openpedcan.urllib.request, "urlopen", urlopen)


def test_fetch_verifies_md5_and_writes_manifest(tmp_path, monkeypatch):
    _fake_bucket(monkeypatch, {"histologies.tsv": b"a\tb\n1\t2\n"})
    got = openpedcan.fetch(tmp_path, files=["histologies.tsv"], release="vX", log=lambda *a: None)
    m = manifest.load(tmp_path / "openpedcan" / "vX")
    assert got["histologies.tsv"].read_bytes() == b"a\tb\n1\t2\n"
    e = m["files"]["histologies.tsv"]
    assert e["sha256"] == hashlib.sha256(b"a\tb\n1\t2\n").hexdigest() and e["retrieved"].endswith("Z")
    assert m["license"] == "CC-BY-4.0" and m["release"] == "vX"


def test_fetch_rejects_md5_mismatch(tmp_path, monkeypatch):
    _fake_bucket(monkeypatch, {"x.tsv": b"data"}, md5_override={"x.tsv": "0" * 32})
    with pytest.raises(ValueError, match="md5 mismatch"):
        openpedcan.fetch(tmp_path, files=["x.tsv"], release="vX", log=lambda *a: None)
    assert not (tmp_path / "openpedcan" / "vX" / "x.tsv").exists()


def test_fetch_skips_existing_valid_file(tmp_path, monkeypatch):
    _fake_bucket(monkeypatch, {"x.tsv": b"data"})
    openpedcan.fetch(tmp_path, files=["x.tsv"], release="vX", log=lambda *a: None)
    seen = []
    openpedcan.fetch(tmp_path, files=["x.tsv"], release="vX", log=seen.append)
    assert seen == ["  have x.tsv"]


# ---------------------------------------------------------------- cohorts

def _hist(rows):
    cols = ["Kids_First_Biospecimen_ID", "sample_id", "cohort", "cancer_group", "composition", "sample_type",
            "experimental_strategy", "RNA_library"]
    return pd.DataFrame(rows, columns=cols)


def test_rna_samples_labels_models_and_drops_adult_and_non_independent():
    h = _hist([
        ["T1", "t1", "TARGET", "Neuroblastoma", "Solid Tissue", "Tumor", "RNA-Seq", "poly-A"],
        ["T2", "t2", "TARGET", "Neuroblastoma", "Solid Tissue", "Tumor", "RNA-Seq", "poly-A"],  # not independent
        ["A1", "a1", "TCGA", "Neuroblastoma", "Solid Tissue", "Tumor", "RNA-Seq", "poly-A"],  # adult cohort
        ["M1", "SK-N-AS", "Maris", "Neuroblastoma", "Derived Cell Line", "Tumor", "RNA-Seq", "poly-A"],
        ["N1", "n1", "GTEx", None, "Solid Tissue", "Normal", "RNA-Seq", "poly-A"],
    ])
    ind = pd.DataFrame({"Kids_First_Biospecimen_ID": ["T1", "A1"]})
    s = cohorts.rna_samples(h, ["T1", "T2", "A1", "M1", "N1"], ind)
    assert dict(zip(s["Kids_First_Biospecimen_ID"], s["kind"])) == {"T1": "tumour", "M1": "model"}


def test_diagnosis_map_reports_unmapped(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text("# c\ncancer_group\toncotree_code\tcanonical_tissue\nNeuroblastoma\tNBL\tpns\nLow-grade glioma\t\tcns\n")
    mapped, unmapped = cohorts.map_groups(["Neuroblastoma", "Low-grade glioma", "Novel thing"],
                                          cohorts.load_diagnosis_map(p))
    assert mapped == {"Neuroblastoma": "NBL"} and unmapped == ["Low-grade glioma", "Novel thing"]


def _concordance_world(seed=0, n_genes=400, n_t=25):
    """Three tumour groups with distinct programs; models = own program + a large
    shared 'culture' offset + noise. One model is mislabelled; one model is in a
    group with no tumours."""
    rng = np.random.default_rng(seed)
    genes = [f"G{i}" for i in range(n_genes)]
    base = rng.normal(5, 1, n_genes)
    prog = {g: rng.normal(0, 1.5, n_genes) for g in ("A", "B", "C")}
    culture = rng.normal(0, 3, n_genes)
    cols, rows = {}, []
    for g in prog:
        for i in range(n_t):
            sid = f"T{g}{i}"
            cols[sid] = base + prog[g] + rng.normal(0, 1, n_genes)
            rows.append([sid, sid, "TARGET", g, "Solid Tissue", "Tumor", "RNA-Seq", "poly-A"])
        for i in range(4):
            sid = f"M{g}{i}"
            cols[sid] = base + prog[g] + culture + rng.normal(0, 1, n_genes)
            rows.append([sid, sid, "Maris", g, "Derived Cell Line", "Tumor", "RNA-Seq", "poly-A"])
    cols["Mbad"] = base + prog["C"] + culture + rng.normal(0, 1, n_genes)  # labelled A, looks like C
    rows.append(["Mbad", "Mbad", "Maris", "A", "Derived Cell Line", "Tumor", "RNA-Seq", "poly-A"])
    cols["Morphan"] = base + culture + rng.normal(0, 1, n_genes)
    rows.append(["Morphan", "Morphan", "Maris", "Z", "Derived Cell Line", "Tumor", "RNA-Seq", "poly-A"])
    expr = pd.DataFrame(cols, index=genes)
    return expr, cohorts.rna_samples(_hist(rows), expr.columns)


def test_concordance_survives_global_culture_offset():
    expr, s = _concordance_world()
    res = cohorts.concordance(expr, s, min_tumours=20, n_var=200, min_tpm_log=0).set_index("biospecimen_id")
    good = res.loc[[f"M{g}{i}" for g in "ABC" for i in range(4)]]
    assert (good["flag"] == "concordant").all() and (good["best_group"] == good["cancer_group"]).all()
    assert res.at["Mbad", "best_group"] == "C" and res.at["Mbad", "flag"] != "concordant"
    assert res.at["Morphan", "flag"] == "no_cohort" and np.isnan(res.at["Morphan", "concordance"])


# ---------------------------------------------------------------- annotations

def _annot_inputs():
    hist = pd.DataFrame({
        "Kids_First_Biospecimen_ID": ["D1", "D2", "D3", "D4", "P1", "R1", "R2"],
        "cohort": ["TARGET"] * 7,
        "cancer_group": ["NB"] * 7,
        "composition": ["Solid Tissue"] * 7,
        "experimental_strategy": ["WGS", "WGS", "WXS", "WGS", "Targeted Sequencing", "RNA-Seq", "RNA-Seq"],
    })
    ind_dna = pd.DataFrame({"Kids_First_Biospecimen_ID": ["D1", "D2", "D3", "D4", "P1"],
                            "experimental_strategy": ["WGS", "WGS", "WXS", "WGS", "Targeted Sequencing"]})
    maf = pd.DataFrame({
        "Hugo_Symbol": ["BRCA2", "BRCA2", "BRCA2", "PARP1"],
        "Variant_Classification": ["Nonsense_Mutation", "Missense_Mutation", "Frame_Shift_Del", "Nonsense_Mutation"],
        "Tumor_Sample_Barcode": ["D1", "D2", "P1", "D3"],  # P1 is a panel sample: excluded
        "HGVSp_Short": ["", "", "", ""],
    })
    cnv = pd.DataFrame({"biospecimen_id": ["D4"], "status": ["deep deletion"], "gene_symbol": ["BRCA2"]})
    fus = pd.DataFrame({"Sample": ["R1"], "Gene1A": ["PARP1"], "Gene1B": ["X"], "Gene2A": [None], "Gene2B": [None]})
    rna = pd.DataFrame({"Kids_First_Biospecimen_ID": ["R1", "R2"], "cancer_group": ["NB", "NB"],
                        "kind": ["tumour", "tumour"]})
    expr = pd.DataFrame({"R1": [1.0, 6.0], "R2": [3.0, 8.0]}, index=["BRCA2", "PARP1"])
    pmtl = pd.DataFrame({"gene_symbol": ["PARP1", "BRCA2"], "pmtl": ["Relevant Molecular Target", None]})
    return hist, rna, expr, ind_dna, maf, cnv, {"D1", "D2", "D3", "D4"}, fus, pmtl


def test_gene_context_counts_lof_deletions_fusions_and_excludes_panels():
    hist, rna, expr, ind_dna, maf, cnv, seen, fus, pmtl = _annot_inputs()
    ctx = annotate.gene_context(["BRCA2", "PARP1"], hist, rna, expr, ind_dna, maf, cnv, seen, fus, pmtl, min_n=2)
    b = ctx.set_index("gene").loc["BRCA2"]
    assert b["n_dna"] == 4 and b["n_lof"] == 1  # nonsense counted; missense and the panel sample not
    assert b["n_deep_del"] == 1 and b["deep_del_fraction"] == 0.25
    p = ctx.set_index("gene").loc["PARP1"]
    assert p["n_fusion"] == 1 and p["fusion_fraction"] == 0.5
    assert p["expr_median"] == 7.0 and p["pmtl"] == "Relevant Molecular Target"


def test_row_summary_reports_baseline_lesion_without_touching_row():
    hist, rna, expr, ind_dna, maf, cnv, seen, fus, pmtl = _annot_inputs()
    ctx = annotate.gene_context(["BRCA2", "PARP1"], hist, rna, expr, ind_dna, maf, cnv, seen, fus, pmtl, min_n=2)
    row = {"lesion_gene": "BRCA2", "drug_target_gene": "PARP1", "drug_name": "OLAPARIB"}
    before = json.dumps(row, sort_keys=True)
    out = annotate.row_summary(row, ctx)
    assert json.dumps(row, sort_keys=True) == before
    assert out["lesion_already_altered"][0]["cancer_group"] == "NB"
    assert out["lesion_max_altered"]["altered_fraction"] == 0.25
    assert out["target_pmtl"] == "Relevant Molecular Target" and "not evidence" in out["note"]
    path = annotate.row_summary(row | {"category": "pathway_suppression"}, ctx)
    assert "lesion_max_altered" not in path and "lesion_note" in path


# ---------------------------------------------------------------- patient side

def test_parquet_cohort_reference_and_tumour_tissue_mismatch(tmp_path, fx, by_id, gene_sets, cfg):
    from npi_pharma.interact.score import score_pair
    from npi_pharma.patient.encode import encode_patient, read_expression

    rng = np.random.default_rng(1)
    genes = list(fx["patient_expr"].index)
    coh = pd.DataFrame(rng.normal(5, 1, (len(genes), 25)), index=genes, columns=[f"S{i}" for i in range(25)])
    coh.to_parquet(tmp_path / "nbl.parquet")
    expr, _ = read_expression(tmp_path / "nbl.parquet")
    p = encode_patient(expr, "S0", "neuroblastoma")
    assert "cohort_relative" in p.flags and "n=24" in p.reference
    rep = score_pair(p, by_id["FIXTURE_LCD_adipose"], by_id["FIXTURE_metformin"], gene_sets, cfg)
    assert "tissue_mismatch" in rep["flags"]
