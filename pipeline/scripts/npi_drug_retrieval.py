"""Which drugs does each NPI resemble? Signature retrieval, the one validated method.

Of everything tried in this project, retrieval is what held up against independent
data: the GSE95640 LCD signature recovered the drug hits an independent published
analysis reported, at AUC 0.75, p = 0.0005 (docs/milestone_GSE95640.md). It answers
a different question from synergy — which is not predictable from signatures
(docs/validation_drugcomb.md) — namely how much an NPI's transcriptional response
resembles a drug's.

  cos > 0   the NPI moves expression the same way the drug does: a MIMIC. Pairing
            them risks duplicating one mechanism rather than adding two.
  cos < 0   the NPI opposes the drug's response.

Read as mechanism only. Similarity is not synergy, and this project measured that
orthogonality does not predict synergy either, so "different, therefore additive"
does not follow.

NPI side: human in-vitro cancer-cell-line signatures only. The muscle, blood and
adipose in-vivo signatures are excluded because the comparison drugs are all
cell-line experiments.

Drug side, in order of preference:
  1. `configs/drug_signatures_geo.yaml` — acute GEO drug-treatment series, which
     cover the platinums, 5-FU and erastin that LINCS Phase II lacks entirely.
  2. LINCS L1000 consensus for the curated drugs it does carry.
  3. A LINCS background panel, used only to place each similarity in a percentile.

    python scripts/npi_drug_retrieval.py --cancer-resource "<...>/Cancer Resource"
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from npi_pharma.store import load_signatures

IN_VITRO_CANCER_CONTEXTS = {"mcf7", "t47d", "mdamb231", "mdamb468", "mcf10a", "u937", "hsc3",
                            "u87", "lovo", "pc3", "lncap"}
# Cell line each NPI signature was measured in, for spotting matched comparisons.
NPI_CELL_LINE = {
    "gse48398_mcf7_heat45c30min": "MCF7", "gse48398_mda231_heat45c30min": "MDA-MB-231",
    "gse48398_mda468_heat45c30min": "MDA-MB-468", "gse48398_mcf10a_heat45c30min": "MCF10A",
    "gse10043_u937_mildhyperthermia41c30min": "U937", "gse75127_hsc3_hyperthermia44c90min": "HSC-3",
    "gse153830_mcf7_bhb10mm": "MCF7", "gse153830_mcf7_glucose_deprivation": "MCF7",
    "gse153830_t47d_bhb25mm": "T47D", "gse153830_t47d_glucose_deprivation": "T47D",
    "gse300765_u87_acidosis_ph64_48h": "U-87", "gse300765_u87_acidosis_ph64_10wk": "U-87",
    "gse300765_u87_hypoxia1pct48h": "U-87", "gse70976_lovo_serumfree96h": "LoVo",
    "gse62673_mcf7_cystine_deprivation24h": "MCF7", "gse62673_mcf7_methionine_deprivation24h": "MCF7",
    "gse62673_pc3_methionine_0um24h": "PC3",
}


def zs(df: pd.DataFrame) -> pd.DataFrame:
    return df.apply(lambda c: (c - c.mean()) / c.std(ddof=0))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cancer-resource", required=True)
    ap.add_argument("--raw-dir", default="data/raw")
    ap.add_argument("--drugs-geo", default="data/processed/signatures/drugs_geo.parquet")
    ap.add_argument("--drugs-lincs", default="data/processed/signatures/lincs_curated.parquet")
    ap.add_argument("--background", default="data/processed/signatures/lincs_background.parquet")
    ap.add_argument("--min-genes", type=int, default=300)
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--out", default="out/npi_drug_retrieval")
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    # NPI side: the human in-vitro cancer-line signatures, plus the GSE62673 arms
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "cnd", Path(__file__).parent / "cancer_npi_drug.py")
    cnd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cnd)
    mods, meta = cnd.modifier_profiles(Path(a.cancer_resource))
    aa = cnd.gse62673_signatures(Path(a.raw_dir))
    keep = [c for c in mods.columns
            if str(meta.loc[c, "context"]).lower() in IN_VITRO_CANCER_CONTEXTS]
    npis = pd.concat([mods[keep], aa], axis=1)
    print(f"NPI signatures (human in-vitro cancer lines): {npis.shape[1]}", file=sys.stderr)

    drugs = {s.sig_id: s for s in load_signatures(a.drugs_geo)}
    src = {k: "GEO acute treatment" for k in drugs}
    for s in load_signatures(a.drugs_lincs):
        drugs[f"{s.sig_id}|LINCS"] = s
        src[f"{s.sig_id}|LINCS"] = "LINCS L1000 consensus"
    D = pd.DataFrame({k: v.as_series() for k, v in drugs.items()})
    B = pd.DataFrame({s.sig_id: s.as_series() for s in load_signatures(a.background)})
    print(f"drug signatures: {D.shape[1]} ({len(load_signatures(a.drugs_geo))} GEO); "
          f"background: {B.shape[1]}", file=sys.stderr)

    # QC first: how much does a drug signature depend on the drug versus the cell
    # line it was measured in? Cross-cell-line retrieval is only interpretable if
    # the same drug agrees between lines.
    gi = D.dropna().index
    Zd = zs(D.loc[gi])
    C = (Zd.T @ Zd) / len(gi)
    C.to_csv(out / "drug_drug_similarity.tsv", sep="\t")
    qc = {"n_shared_genes": int(len(gi))}
    pairs = [("erastin|LNCaP", "erastin|PC3", "same drug, two cell lines"),
             ("oxaliplatin|HCT116", "5-fluorouracil|HCT116", "same cell line, two DNA-damaging drugs"),
             ("sirolimus|LINCS", "everolimus|LINCS", "two rapalogs, LINCS consensus"),
             ("cisplatin|MDA-MB-231", "oxaliplatin|HCT116", "two platinums, different cell lines")]
    for x, y, why in pairs:
        if x in C.index and y in C.columns:
            qc[why] = float(C.at[x, y])
    print("\nQC, drug-signature similarity:", file=sys.stderr)
    for k, v in qc.items():
        if k != "n_shared_genes":
            print(f"  {k}: {v:+.3f}", file=sys.stderr)

    rows = []
    for npi in npis.columns:
        n = npis[npi].dropna()
        for target, label in ((D, "drug"), (B, "background")):
            g = n.index.intersection(target.index)
            if len(g) < a.min_genes:
                continue
            x = (n[g] - n[g].mean()) / n[g].std(ddof=0)
            Z = zs(target.loc[g])
            cos = (Z.T.fillna(0) @ x.fillna(0)) / len(g)
            for k, v in cos.items():
                rows.append({"npi": npi, "npi_cell_line": NPI_CELL_LINE.get(npi),
                             "target": k, "kind": label, "similarity": float(v),
                             "n_genes": len(g)})
    t = pd.DataFrame(rows)
    bg = t[t["kind"] == "background"]
    res = t[t["kind"] == "drug"].copy()
    # percentile of each similarity against that NPI's own background distribution
    res["background_percentile"] = [
        float((bg.loc[bg["npi"] == r.npi, "similarity"] < r.similarity).mean())
        for r in res.itertuples()]
    res["drug_source"] = res["target"].map(src)
    res["drug_cell_line"] = res["target"].map(
        {k: (v.meta or {}).get("cell_line") for k, v in drugs.items()})
    res["same_cell_line"] = [
        bool(r.npi_cell_line and r.drug_cell_line and r.npi_cell_line == r.drug_cell_line)
        for r in res.itertuples()]
    res = res.sort_values(["npi", "similarity"], ascending=[True, False])
    res.to_csv(out / "npi_drug_similarity.tsv", sep="\t", index=False)

    matched = res[res["same_cell_line"]]
    summary = {
        "n_npis": int(res["npi"].nunique()), "n_drug_signatures": int(res["target"].nunique()),
        "n_background": int(B.shape[1]),
        "interpretation": "cos > 0: NPI mimics the drug (mechanism duplication risk). "
                          "cos < 0: opposes. Similarity is not synergy.",
        "same_cell_line_comparisons": matched[["npi", "target", "similarity", "background_percentile"]]
            .to_dict(orient="records"),
        "drug_signature_qc": qc,
        "qc_warning": "A drug signature depends more on the cell line it was measured in than on the "
                      "drug: same cell line with two different DNA-damaging drugs reaches ~0.71 while "
                      "the same drug in two cell lines reaches only ~0.18. Cross-cell-line similarity "
                      "is therefore weakly interpretable, and the LINCS consensus signatures (medians "
                      "over many cell lines) are the ones the validated GSE95640 retrieval used.",
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float))

    print(f"\n{res['npi'].nunique()} NPIs x {res['target'].nunique()} drug signatures "
          f"({B.shape[1]}-compound background)\n")
    show = ["target", "similarity", "background_percentile", "drug_source", "n_genes"]
    for npi, g in res.groupby("npi"):
        print(f"== {npi}  [{NPI_CELL_LINE.get(npi, '?')}]")
        print(g.head(a.top)[show].to_string(index=False))
        print("   most opposed: " + ", ".join(
            f"{r.target} {r.similarity:+.3f}" for r in g.tail(2).itertuples()))
    if len(matched):
        print("\n== same-cell-line comparisons (NPI and drug measured in the same line)")
        print(matched[["npi", "npi_cell_line", "target", "similarity", "background_percentile"]]
              .to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
