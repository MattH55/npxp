"""NPI x drug similarity across every drug signature available, with reliability flags.

Consolidates the drug panels this project has built and ranks each NPI against all
of them. Weak signatures are **included and labelled**, not withheld, so a result
can be read with its confidence attached rather than silently dropped.

Reliability, from this project's own measurements:

  high      LINCS consensus — a median over a median of 7 cell lines per compound.
  medium    GEO consensus with cross-series agreement >= 0.30. Validation against
            LINCS on 10 drugs put agreement at 0.39 median, and cross-series
            agreement predicts that at r = 0.70 (docs/drug_consensus.md).
  low       GEO consensus with cross-series agreement 0.10-0.30.
  very low  GEO consensus below 0.10, or a single-series signature. Single-cell-line
            signatures are dominated by the cell line rather than the drug: the same
            drug in two lines agrees at 0.18 while two different drugs in one line
            agree at 0.71 (docs/npi_drug_retrieval.md).

What the score means: cos > 0 says the NPI moves expression the way the drug does
(a mimic, so pairing them risks duplicating one mechanism); cos < 0 says it opposes.
Similarity is NOT synergy — this project measured that neither similarity nor
orthogonality predicts it (docs/validation_drugcomb.md).

    python scripts/rank_npi_drug_all.py --cancer-resource "<...>/Cancer Resource"
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from npi_drug_retrieval import IN_VITRO_CANCER_CONTEXTS, NPI_CELL_LINE, zs  # noqa: E402

from npi_pharma.store import load_signatures  # noqa: E402

PANELS = [
    ("data/processed/signatures/lincs_all.parquet", "LINCS consensus", "high"),
    ("data/processed/signatures/drugs_consensus.parquet", "GEO consensus", None),
    ("data/processed/signatures/drugs_geo.parquet", "GEO single series", "very low"),
]


def independent_accessions(accessions: list[str], adjacent_within: int = 5) -> bool:
    """Do these GEO accessions represent genuinely independent studies?

    Two arms of one accession share protocol, batch and operator, so they are not
    independent however well they agree. Consecutive accession numbers are almost
    always companion series of a single submission (GSE59296 and GSE59297 are one
    topotecan study split by platform), so they are treated the same way.
    """
    uniq = sorted({a for a in accessions if a})
    if len(uniq) < 2:
        return False
    nums = sorted(int(m.group(1)) for a in uniq if (m := re.match(r"GSE(\d+)$", a)))
    if len(nums) == len(uniq) and all(b - a <= adjacent_within for a, b in zip(nums, nums[1:])):
        return False
    return True


def reliability_for(drug: str, qc: dict, accessions: list[str] | None = None
                    ) -> tuple[str, float | None, int | None]:
    """Label a GEO consensus from its cross-series agreement.

    A consensus whose arms all come from ONE accession is not independent evidence
    however well they agree -- two cell lines of the same study share protocol,
    batch and operator -- so it is capped at "very low" regardless of its number.
    """
    entry = qc.get(drug) or qc.get(drug.replace("-", " ")) or {}
    r = entry.get("median_cross_series_cosine")
    n = entry.get("n_series")
    if accessions is not None and not independent_accessions(accessions):
        return "very low", (float(r) if r is not None else None), n
    if r is None:
        return "very low", None, n
    if r >= 0.30:
        return "medium", float(r), n
    if r >= 0.10:
        return "low", float(r), n
    return "very low", float(r), n


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cancer-resource", required=True)
    ap.add_argument("--raw-dir", default="data/raw")
    ap.add_argument("--consensus-qc", default="out/drug_consensus/qc.json")
    ap.add_argument("--min-genes", type=int, default=300)
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--out", default="out/npi_drug_all")
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    qc = json.loads(Path(a.consensus_qc).read_text()) if Path(a.consensus_qc).exists() else {}

    import importlib.util

    spec = importlib.util.spec_from_file_location("cnd", Path(__file__).parent / "cancer_npi_drug.py")
    cnd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cnd)
    mods, meta = cnd.modifier_profiles(Path(a.cancer_resource))
    aa = cnd.gse62673_signatures(Path(a.raw_dir))
    keep = [c for c in mods.columns if str(meta.loc[c, "context"]).lower() in IN_VITRO_CANCER_CONTEXTS]
    npis = pd.concat([mods[keep], aa], axis=1)

    frames, info = {}, {}
    for path, source, fixed in PANELS:
        p = Path(path)
        if not p.exists():
            print(f"  (missing: {path})", file=sys.stderr)
            continue
        for s in load_signatures(p):
            drug = (s.meta or {}).get("drug_id", s.sig_id)
            key = f"{s.sig_id}|{source}"
            frames[key] = s.as_series()
            rel, cross, nser = ((fixed, None, None) if fixed else
                                reliability_for(drug, qc, list(s.source_accessions or [])))
            info[key] = {"drug": drug, "source": source, "reliability": rel,
                         "cross_series_agreement": cross, "n_series": nser,
                         "cell_line": (s.meta or {}).get("cell_line")}
    D = pd.DataFrame(frames)
    print(f"{npis.shape[1]} NPI signatures x {D.shape[1]} drug signatures", file=sys.stderr)
    counts = pd.Series([v["reliability"] for v in info.values()]).value_counts()
    print("drug signatures by reliability: " + ", ".join(f"{k} {v}" for k, v in counts.items()),
          file=sys.stderr)

    rows = []
    for npi in npis.columns:
        n = npis[npi].dropna()
        g = n.index.intersection(D.index)
        if len(g) < a.min_genes:
            continue
        x = (n[g] - n[g].mean()) / n[g].std(ddof=0)
        Z = zs(D.loc[g])
        cos = (Z.T.fillna(0) @ x.fillna(0)) / len(g)
        for k, v in cos.items():
            rows.append({"npi": npi, "npi_cell_line": NPI_CELL_LINE.get(npi),
                         "drug_signature": k, "similarity": float(v), "n_genes": len(g), **info[k]})
    t = pd.DataFrame(rows)
    # percentile within each NPI, against the high-reliability LINCS panel only
    ref = t[t["reliability"] == "high"]
    t["percentile_vs_lincs_panel"] = [
        float((ref.loc[ref["npi"] == r.npi, "similarity"] < r.similarity).mean())
        for r in t.itertuples()]
    t = t.sort_values(["npi", "similarity"], ascending=[True, False])
    t.to_csv(out / "npi_drug_all.tsv", sep="\t", index=False)

    summary = {
        "n_npis": int(t["npi"].nunique()), "n_drug_signatures": int(t["drug_signature"].nunique()),
        "by_reliability": {k: int(v) for k, v in counts.items()},
        "interpretation": "cos > 0: the NPI mimics the drug (mechanism duplication risk). "
                          "Similarity is not synergy.",
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float))

    show = ["drug_signature", "similarity", "percentile_vs_lincs_panel", "reliability"]
    print(f"\n{t['npi'].nunique()} NPIs scored against {t['drug_signature'].nunique()} "
          f"drug signatures\n")
    for npi, g in t.groupby("npi"):
        print(f"== {npi}  [{NPI_CELL_LINE.get(npi, '?')}]")
        print(g.head(a.top)[show].to_string(index=False))
    # and the curated drugs specifically, whatever their rank
    curated = t[t["drug"].isin(["cisplatin", "carboplatin", "oxaliplatin", "5-fluorouracil",
                                "mitomycin-c", "doxorubicin", "paclitaxel", "temozolomide",
                                "metformin", "erastin"])]
    curated.to_csv(out / "curated_drugs.tsv", sep="\t", index=False)
    print(f"\ncurated-drug rows written to {out}/curated_drugs.tsv ({len(curated)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
