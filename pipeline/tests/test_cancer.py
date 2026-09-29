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


def test_probe_map_from_platform_table(monkeypatch, tmp_path):
    """Both annotation shapes GEO serves: a symbol column, and Affymetrix gene_assignment."""
    from npi_pharma.ingest import fetch

    symbol_table = (
        "^PLATFORM = GPL1\n!platform_table_begin\n"
        "ID\tCONTROL_TYPE\tGENE_SYMBOL\tGENE_NAME\n"
        "p1\t\tMTOR\tmechanistic target\n"
        "p2\tpos\t\t\n"                       # no symbol -> dropped
        "p3\t\tRPTOR /// RPTOR2\tx\n"          # first of /// kept
    )
    assign_table = (
        "ID\tprobeset_id\tgene_assignment\n"
        "q1\tq1\tNM_001005484 // SAMD11 // sterile alpha motif // 1p36\n"
        "q2\tq2\t---\n"
        "q3\tq3\tNM_000546 // TP53 // tumor protein p53 /// NM_0001 // TP53B // other\n"
    )
    pages = {"A": symbol_table, "B": assign_table}
    monkeypatch.setattr(fetch, "_get", lambda url, timeout=180: pages[url[-1]].encode())
    monkeypatch.setattr(fetch, "PLATFORM_TABLE_URL", "http://x/{gpl}A")
    out = fetch.probe_map_from_platform_table("GPL1", tmp_path / "a.tsv")
    assert out.read_text().splitlines() == ["p1\tMTOR", "p3\tRPTOR"]

    monkeypatch.setattr(fetch, "PLATFORM_TABLE_URL", "http://x/{gpl}B")
    out = fetch.probe_map_from_platform_table("GPL2", tmp_path / "b.tsv")
    assert out.read_text().splitlines() == ["q1\tSAMD11", "q3\tTP53"]

    monkeypatch.setattr(fetch, "PLATFORM_TABLE_URL", "http://x/{gpl}C")
    pages["C"] = "ID\tRANGE_START\tGB_ACC\nr1\t1\tNR_046018\n"
    with pytest.raises(ValueError, match="no gene-symbol column"):
        fetch.probe_map_from_platform_table("GPL3", tmp_path / "c.tsv")


def test_consensus_arm_detection_rejects_resistance_and_combos():
    import importlib.util
    from pathlib import Path

    spec = importlib.util.spec_from_file_location(
        "bdc", Path(__file__).parents[1] / "scripts" / "build_drug_consensus.py")
    bdc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bdc)

    idx = [f"GSM{i}" for i in range(6)]
    # a clean acute design is detected
    ok = pd.DataFrame({"treatment": ["DMSO", "DMSO", "DMSO", "cisplatin 10 uM",
                                     "cisplatin 10 uM", "cisplatin 10 uM"]}, index=idx)
    arms = bdc.detect_arms(ok, "cisplatin")
    assert arms is not None and len(arms[0]) == 3 and len(arms[1]) == 3 and arms[2] == "treatment"

    # a resistant-vs-parental comparison is not a drug response
    res = pd.DataFrame({"resistance": ["parental", "parental", "parental",
                                       "cisplatin-resistant", "cisplatin-resistant",
                                       "cisplatin-resistant"]}, index=idx)
    assert bdc.detect_arms(res, "cisplatin") is None

    # a combination arm must not be taken for the drug alone
    combo = pd.DataFrame({"treatment": ["control", "control", "control",
                                        "cisplatin + olaparib", "cisplatin + olaparib",
                                        "cisplatin + olaparib"]}, index=idx)
    assert bdc.detect_arms(combo, "cisplatin") is None

    # too few replicates
    small = pd.DataFrame({"treatment": ["control", "cisplatin"]}, index=idx[:2])
    assert bdc.detect_arms(small, "cisplatin") is None


def _bdc():
    import importlib.util
    from pathlib import Path

    spec = importlib.util.spec_from_file_location(
        "bdc2", Path(__file__).parents[1] / "scripts" / "build_drug_consensus.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_match_columns_uses_geo_labels_and_refuses_ambiguity():
    bdc = _bdc()
    samples = pd.DataFrame(
        {"title": ["ctrl rep 1", "ctrl rep 2", "drug rep 1"],
         "description": ["Library name: NC2_1", "Library name: NC2_2", "Library name: Bo2_1"],
         "supplementary_file_1": ["NONE", "NONE", "NONE"]},
        index=["GSM1", "GSM2", "GSM3"])
    cols = ["gene", "NC2_1", "NC2_2", "Bo2_1"]
    got = bdc.match_columns(samples, ["GSM1", "GSM3"], cols)
    assert got == {"GSM1": "NC2_1", "GSM3": "Bo2_1"}

    # a label embedded in a longer column name still resolves
    cols2 = ["gene_name", "run_NC2_1_count", "run_NC2_2_count", "run_Bo2_1_count"]
    assert bdc.match_columns(samples, ["GSM1"], cols2) == {"GSM1": "run_NC2_1_count"}

    # nothing to match -> refuse rather than guess
    assert bdc.match_columns(samples, ["GSM1"], ["gene", "sampleA", "sampleB"]) is None

    # an ambiguous label (matches two columns) -> refuse
    amb = pd.DataFrame({"title": ["rep1"], "description": ["Library name: A"],
                        "supplementary_file_1": ["NONE"]}, index=["GSM9"])
    assert bdc.match_columns(amb, ["GSM9"], ["A_1", "A_2"]) is None


def test_labels_for_strips_library_name_prefix():
    bdc = _bdc()
    samples = pd.DataFrame(
        {"title": ["MDA-MB-231 cisplatin rep 1"], "description": ["Library name: Bo2_1"],
         "supplementary_file_1": ["ftp://ftp.ncbi.nlm.nih.gov/x/GSM1_counts.txt.gz"]},
        index=["GSM1"])
    labs = bdc._labels_for(samples, "GSM1")
    assert "Bo2_1" in labs                      # prefix stripped
    assert "Library name: Bo2_1" in labs        # and the raw field kept
    assert "GSM1_counts" in labs                # supplementary basename, extensions removed


def test_independent_accessions_rejects_one_study_and_companion_series():
    import importlib.util
    from pathlib import Path

    spec = importlib.util.spec_from_file_location(
        "rnda", Path(__file__).parents[1] / "scripts" / "rank_npi_drug_all.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)

    assert m.independent_accessions(["GSE121689", "GSE146354"]) is True   # far apart: two studies
    assert m.independent_accessions(["GSE232034"]) is False               # one accession
    assert m.independent_accessions(["GSE232034", "GSE232034"]) is False  # same, twice
    assert m.independent_accessions(["GSE59296", "GSE59297"]) is False    # companion series
    assert m.independent_accessions([]) is False
    # three consecutive is still one submission; a distant third makes it independent
    assert m.independent_accessions(["GSE100", "GSE101", "GSE102"]) is False
    assert m.independent_accessions(["GSE100", "GSE101", "GSE90000"]) is True


def test_reliability_bands_follow_cross_series_agreement():
    import importlib.util
    from pathlib import Path

    spec = importlib.util.spec_from_file_location(
        "rnda2", Path(__file__).parents[1] / "scripts" / "rank_npi_drug_all.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)

    qc = {"good": {"median_cross_series_cosine": 0.45, "n_series": 4},
          "mid": {"median_cross_series_cosine": 0.20, "n_series": 3},
          "poor": {"median_cross_series_cosine": 0.02, "n_series": 2}}
    ind = ["GSE1", "GSE9999"]
    assert m.reliability_for("good", qc, ind)[0] == "medium"
    assert m.reliability_for("mid", qc, ind)[0] == "low"
    assert m.reliability_for("poor", qc, ind)[0] == "very low"
    assert m.reliability_for("unknown", qc, ind)[0] == "very low"
    # a strong number from a single study is still capped
    assert m.reliability_for("good", qc, ["GSE1", "GSE1"])[0] == "very low"


def test_match_columns_by_assignment_solves_arm_labelled_tables():
    """Count tables whose columns are arm labels, not GSM names, still map.

    Each case below is a real GEO layout: the assignment has to survive
    abbreviation (Veh/Vehicle), run-together labels (SNB19ACT), a replicate
    number written as _R1, and a control arm identified only by what is left over.
    """
    bdc = _bdc()

    def samples(titles):
        return pd.DataFrame({"title": titles, "description": [""] * len(titles),
                             "supplementary_file_1": ["NONE"] * len(titles)},
                            index=[f"GSM{i}" for i in range(len(titles))])

    s = samples(["CAMA1 DMSO RNAseq replicate 1", "CAMA1 abemaciclib RNAseq replicate 1",
                 "CAMA1 fulvestrant RNAseq replicate 1"])
    assert bdc.match_columns_by_assignment(
        s, list(s.index), ["CAMA-1_DMSO_1", "CAMA-1_Abema_1", "CAMA-1_Fulv_1"]) == {
            "GSM0": "CAMA-1_DMSO_1", "GSM1": "CAMA-1_Abema_1", "GSM2": "CAMA-1_Fulv_1"}

    # abbreviated labels, and "_R1" against a spelled-out "Replicate 1"
    s = samples(["Vehicle R", "Palbociclib R", "DMSO Replicate 1", "VTP_WM119 Replicate 1"])
    assert bdc.match_columns_by_assignment(
        s, list(s.index), ["Veh_R", "Palbo_R", "DMSO_R1", "VTP_WM119_R1"]) == {
            "GSM0": "Veh_R", "GSM1": "Palbo_R", "GSM2": "DMSO_R1", "GSM3": "VTP_WM119_R1"}

    # the control column carries no word for "control": it is the one left over
    s = samples(["SNB19 cells,  Control", "SNB19 cells,  ACT001", "SNB19 cells,  Stattic"])
    assert bdc.match_columns_by_assignment(
        s, list(s.index), ["SNB19", "SNB19ACT", "SNB19Sta"]) == {
            "GSM0": "SNB19", "GSM1": "SNB19ACT", "GSM2": "SNB19Sta"}

    # counts and TPM of the same samples are two value types, not two labels
    s = samples(["DMSO_1", "DAL_1"])
    assert bdc.match_columns_by_assignment(
        s, list(s.index), ["DMSO_1_tpm", "DAL_1_tpm", "DMSO_1_count", "DAL_1_count"]) == {
            "GSM0": "DMSO_1_count", "GSM1": "DAL_1_count"}


def test_match_columns_by_assignment_refuses_when_not_decisive():
    bdc = _bdc()

    def samples(titles):
        return pd.DataFrame({"title": titles, "description": [""] * len(titles),
                             "supplementary_file_1": ["NONE"] * len(titles)},
                            index=[f"GSM{i}" for i in range(len(titles))])

    # two samples that cannot be told apart: swapping them scores the same
    s = samples(["treated", "treated"])
    assert bdc.match_columns_by_assignment(s, list(s.index), ["T_a", "T_b"]) is None

    # labels that name nothing in the columns
    s = samples(["alpha", "beta"])
    assert bdc.match_columns_by_assignment(s, list(s.index), ["colX", "colY"]) is None

    # one sample gives the joint constraint nothing to work with
    s = samples(["DMSO Replicate 1", "other"])
    assert bdc.match_columns_by_assignment(s, ["GSM0"], ["DMSO_R1", "other"]) is None

    # fewer columns than samples
    s = samples(["DMSO Replicate 1", "VTP Replicate 1"])
    assert bdc.match_columns_by_assignment(s, list(s.index), ["DMSO_R1"]) is None


def test_consensus_driver_slug_makes_safe_part_names():
    import importlib.util
    from pathlib import Path

    spec = importlib.util.spec_from_file_location(
        "bdca", Path(__file__).parents[1] / "scripts" / "build_drug_consensus_all.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)

    assert m.slug("5-fluorouracil") == "5_fluorouracil"
    assert m.slug("mitomycin C") == "mitomycin_c"
    assert m.slug("actinomycin D") == "actinomycin_d"
    # distinct drugs must not collide on one part directory
    names = [m.slug(d) for d in m.MISSING_FROM_LINCS]
    assert len(set(names)) == len(names)
