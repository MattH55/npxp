import numpy as np
import pandas as pd
import pytest

from npi_pharma.cancer.programs import center_within, drug_specific, resistance_programs
from npi_pharma.cancer.sensitize import percentile, score_matrix


def _world(seed=0, n=300, g=300):
    rng = np.random.default_rng(seed)
    lines = [f"ACH-{i:06d}" for i in range(n)]
    genes = [f"G{i}" for i in range(g)]
    lineage = pd.Series(np.repeat(["lung", "breast", "skin"], n // 3), index=lines)
    x = pd.DataFrame(rng.normal(size=(n, g)), index=lines, columns=genes)
    x.loc[lineage == "lung"] += 2.0  # lineage shift on every gene
    auc = pd.DataFrame({
        "drugA": 0.8 * x["G0"] + rng.normal(0, 0.5, n),         # G0 high -> resistant
        "drugB": -0.8 * x["G1"] + rng.normal(0, 0.5, n),        # G1 high -> sensitive
        "tissueonly": (lineage == "lung").astype(float) * 3 + rng.normal(0, 0.1, n),
    }, index=lines)
    return x, auc, lineage


def test_center_within_drops_small_groups():
    df = pd.DataFrame({"v": [1.0, 3.0, 5.0, 10.0]}, index=list("abcd"))
    g = pd.Series(["x", "x", "x", "y"], index=list("abcd"))
    c = center_within(df, g, min_group=2)
    assert list(c.index) == ["a", "b", "c"] and c["v"].tolist() == [-2.0, 0.0, 2.0]


def test_resistance_programs_recover_planted_genes_and_ignore_lineage():
    x, auc, lineage = _world()
    progs, n = resistance_programs(x, auc, lineage, min_lines=50)
    assert progs["drugA"].idxmax() == "G0" and progs["drugB"].idxmin() == "G1"
    # AUC driven only by lineage leaves no program after lineage centring
    assert progs["tissueonly"].abs().max() < 0.25
    assert n["drugA"] == 300


def test_drug_specific_removes_shared_axis():
    genes = [f"G{i}" for i in range(50)]
    rng = np.random.default_rng(1)
    shared = rng.normal(size=50)
    progs = pd.DataFrame({f"d{i}": shared + 0.3 * rng.normal(size=50) for i in range(20)}, index=genes)
    spec = drug_specific(progs)
    bg = progs.mean(axis=1).to_numpy()
    assert all(abs(spec[c].to_numpy() @ bg) < 1e-8 for c in spec)  # orthogonal to the shared axis
    assert np.linalg.norm(spec.mean(axis=1)) < 1e-8


def test_score_sign_and_percentile():
    genes = [f"G{i}" for i in range(300)]
    rng = np.random.default_rng(2)
    r = pd.DataFrame({"resist_up": rng.normal(size=300), "other": rng.normal(size=300)}, index=genes)
    npi = pd.DataFrame({"moves_to_sensitive": -r["resist_up"] + 0.1 * rng.normal(size=300)}, index=genes)
    s = score_matrix(npi, r)
    assert s.loc["moves_to_sensitive", "resist_up"] > 0.9
    assert percentile(s).loc["moves_to_sensitive", "resist_up"] == 1.0


def test_residual_spearman_removes_main_effects():
    import importlib.util
    from pathlib import Path

    spec = importlib.util.spec_from_file_location(
        "vm", Path(__file__).parents[1] / "scripts" / "validate_monotherapy.py")
    vm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(vm)
    rng = np.random.default_rng(0)
    drugs, cells = [f"d{i}" for i in range(20)], [f"c{i}" for i in range(20)]
    de = dict(zip(drugs, rng.normal(size=20)))
    ce = dict(zip(cells, rng.normal(size=20)))

    rows = []
    for d in drugs:
        for c in cells:
            pair = rng.normal()              # genuinely pair-specific, not additive
            rows.append({"drug": d, "cell": c, "additive": de[d] + ce[c],
                         "pair_only": pair, "auc": de[d] + ce[c] + pair})
    t = pd.DataFrame(rows)

    # a purely additive outcome leaves a numerically zero residual: no signal to find
    add = t.assign(auc=t["additive"])
    v = add["auc"]
    resid = v - add.groupby("drug")["auc"].transform("mean") - add.groupby("cell")["auc"].transform("mean") + v.mean()
    assert resid.abs().max() < 1e-12

    # a feature that knows only the pair-specific term scores ~0 raw but ~1 on the residual
    raw, res, n = vm.residual_spearman(t, "pair_only", "auc", "drug", "cell")
    assert abs(raw) < 0.75 and res > 0.95 and n == 400
    # a feature that knows only the main effects scores high raw but nothing on the residual
    raw2, res2, _ = vm.residual_spearman(t, "additive", "auc", "drug", "cell")
    assert raw2 > 0.5 and abs(res2) < 0.1


def test_bliss_expected_properties():
    from npi_pharma.cancer.efficacy import auc_to_inhibition, bliss_expected, excess_over_bliss

    assert bliss_expected(0.3, 0.5) == pytest.approx(0.65)
    assert bliss_expected(0.0, 0.5) == pytest.approx(0.5)          # an inert agent adds nothing
    assert bliss_expected(-0.2, 0.5) == pytest.approx(0.5)         # growth-promoting is clipped, not credited
    assert bliss_expected(1.0, 0.4) == pytest.approx(1.0)          # saturates
    assert bliss_expected(0.5, 0.3) == pytest.approx(bliss_expected(0.3, 0.5))  # symmetric
    # monotone in each argument
    assert bliss_expected(0.4, 0.5) > bliss_expected(0.3, 0.5)
    np.testing.assert_allclose(bliss_expected(np.array([0.0, 0.5]), np.array([0.5, 0.5])), [0.5, 0.75])
    # excess is measured minus expected; zero when the measurement is exactly additive
    assert excess_over_bliss(0.65, 0.3, 0.5) == pytest.approx(0.0)
    assert excess_over_bliss(0.8, 0.3, 0.5) == pytest.approx(0.15)
    assert auc_to_inhibition(0.7) == pytest.approx(0.3)
    assert auc_to_inhibition(1.2) == pytest.approx(0.0)            # AUC above 1 clipped to no effect


def test_rank_combinations_drops_cell_lines_without_npi_measurement():
    from npi_pharma.cancer.efficacy import rank_combinations

    npi = pd.Series({"RKO": 0.8, "HCT116": 0.0})
    drugs = pd.DataFrame({"a": [0.5, 0.5], "b": [0.1, np.nan]}, index=["RKO", "OTHER"])
    out = rank_combinations(npi, drugs)
    assert set(out["cell_line_id"]) == {"RKO"}          # HCT116 has no drug data, OTHER has no NPI value
    assert set(out["drug"]) == {"a", "b"}               # NaN drug values are dropped, not imputed
    assert out.iloc[0]["drug"] == "a"
    assert out.iloc[0]["expected_combined_inhibition"] == pytest.approx(0.9)


def test_drug_geo_gene_index_from_compound_id():
    from npi_pharma.ingest.drug_geo import _gene_index

    df = pd.DataFrame({"test_id": ["ENSG00000000003_TSPAN6", "ENSG00000000005_TNMD"], "NT1": [1, 2]})
    g = _gene_index(df, {"id_col": "test_id", "gene_from_id": "_"})
    assert list(g) == ["TSPAN6", "TNMD"]
    g2 = _gene_index(pd.DataFrame({"gene_name": ["A", "B"]}), {"gene_col": "gene_name"})
    assert list(g2) == ["A", "B"]
    with pytest.raises(ValueError, match="cannot find gene symbols"):
        _gene_index(pd.DataFrame({"x": [1]}), {})


def test_drug_geo_build_filters_and_labels(tmp_path):
    from npi_pharma.ingest.drug_geo import PROV_GEO_DRUG, build_from_record
    from npi_pharma.model import DRUG

    rng = np.random.default_rng(0)
    n = 3000
    genes = [f"G{i}" for i in range(n)]
    base = rng.integers(200, 2000, size=(n, 1)) * np.ones((1, 6))
    base = base + rng.integers(0, 40, size=(n, 6))
    base[:20, 3:] *= 5          # planted up in treated
    base[20:40, 3:] = base[20:40, 3:] // 5   # planted down
    df = pd.DataFrame(base, columns=["c1", "c2", "c3", "t1", "t2", "t3"])
    df.insert(0, "gene_name", genes)
    df.loc[n - 1, "gene_name"] = ""          # blank symbol must be dropped
    path = tmp_path / "expr.tsv"
    df.to_csv(path, sep="\t", index=False)
    rec = {"drug_id": "testdrug", "accession": "GSE1", "cell_line": "XYZ", "cell_line_id": None,
           "curated_cell_line": False, "condition": "1 uM, 24 h", "file": "expr.tsv",
           "format": {"sep": "\t", "gene_col": "gene_name"},
           "control": ["c1", "c2", "c3"], "treated": ["t1", "t2", "t3"]}
    sig = build_from_record(rec, tmp_path, min_genes=100)
    s = sig.as_series()
    assert sig.kind == DRUG and sig.provenance == PROV_GEO_DRUG and sig.sig_id == "testdrug|XYZ"
    assert sig.meta["n_control"] == 3 and sig.meta["drug_id"] == "testdrug"
    assert "" not in s.index
    assert set(s.nlargest(10).index) <= {f"G{i}" for i in range(20)}
    assert set(s.nsmallest(10).index) <= {f"G{i}" for i in range(20, 40)}
    rec_bad = dict(rec, control=["nope"])
    with pytest.raises(ValueError, match="columns not in file"):
        build_from_record(rec_bad, tmp_path, min_genes=100)
