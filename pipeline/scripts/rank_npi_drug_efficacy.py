"""Rank NPI x drug per cell line by expected combination EFFICACY.

This applies the literature's finding rather than working around it: each agent's
own response in the cell line is the strongest available predictor of what a
combination does, and total kill is predictable where excess-over-additivity is
not (docs/efficacy_from_monotherapy.md, docs/validation_drugcomb.md).

  NPI side   configs/npi_monotherapy.yaml -- curated published inhibition of the
             NPI alone in a named cell line. Absent entries are reported as
             missing coverage and drop that cell line; nothing is imputed.
  drug side  PRISM 19Q4 and/or GDSC fitted AUC for that cell line, rescaled to
             inhibition by 1 - AUC.
  score      Bliss independence, 1 - (1-npi)(1-drug), which on DrugComb ranks
             measured combination efficacy at within-cell-line Spearman 0.752
             with no fitting.

The output is a ranking of expected total kill. It is NOT a synergy prediction:
measured efficacy and Bliss synergy correlate 0.060 in DrugComb.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from npi_pharma.cancer.efficacy import auc_to_inhibition, bliss_expected
from npi_pharma.cancer.programs import load_gdsc_auc, load_prism_auc


def load_npi_monotherapy(path: str | Path) -> pd.DataFrame:
    d = yaml.safe_load(Path(path).read_text()) or {}
    rows = []
    for r in d.get("npi_monotherapy", []):
        rows.append({"npi_id": r["npi_id"], "cell_line_id": r["cell_line_id"], "assay": r.get("assay"),
                     "metric": r.get("metric"), "inhibition": float(r["inhibition"]),
                     "temperature_C": (r.get("condition") or {}).get("temperature_C"),
                     "source": (r.get("source") or {}).get("pmid")})
    return pd.DataFrame(rows)


def to_depmap_id(cell_line_id: str, models: pd.DataFrame) -> str | None:
    """Map a curated id to a DepMap ModelID.

    Accepts a ModelID as-is, or this project's CCLE-style ids where the name is
    followed by the tissue (RKO_LARGE_INTESTINE, MDAMB231_BREAST): the part before
    the first underscore is matched against DepMap's stripped cell-line names.
    """
    import re

    norm = lambda s: re.sub(r"[^A-Z0-9]", "", str(s).upper())
    if str(cell_line_id).startswith("ACH-"):
        return str(cell_line_id)
    table = {norm(r.StrippedCellLineName): r.ModelID for r in models.itertuples()
             if isinstance(r.StrippedCellLineName, str)}
    key = norm(str(cell_line_id).split("_")[0])
    return table.get(key)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--npi-monotherapy", default="configs/npi_monotherapy.yaml")
    ap.add_argument("--drug-sets", default="configs/drug_sets.yaml")
    ap.add_argument("--drug-set", help="restrict to a named drug set")
    ap.add_argument("--depmap-dir", default="data/raw/depmap")
    ap.add_argument("--assay", default="clonogenic_survival_10d",
                    help="which NPI assay to prefer; falls back to any assay for that cell line")
    ap.add_argument("--out", default="out/npi_drug_efficacy")
    a = ap.parse_args(argv)
    out, dm = Path(a.out), Path(a.depmap_dir)
    out.mkdir(parents=True, exist_ok=True)

    npi = load_npi_monotherapy(a.npi_monotherapy)
    if npi.empty:
        print("no curated NPI monotherapy entries: nothing can be ranked", file=sys.stderr)
        return 1
    models = pd.read_csv(dm / "Model.csv")
    name_of = models.set_index("ModelID")["StrippedCellLineName"]
    npi["model_id"] = npi["cell_line_id"].map(lambda c: to_depmap_id(c, models))
    unmapped = sorted(npi.loc[npi["model_id"].isna(), "cell_line_id"].unique())
    if unmapped:
        print(f"  no DepMap cell line for: {', '.join(unmapped)}", file=sys.stderr)

    sources = {"PRISM": load_prism_auc(dm / "prism19q4_secondary_dose_response.csv")}
    for p in sorted(dm.glob("GDSC*_fitted_dose_response*.csv")):
        sources[p.name.split("_")[0]] = load_gdsc_auc(p, models)
    # map PRISM's depmap_id index and GDSC's ModelID index onto DepMap ids (both already are)
    drug_sets = yaml.safe_load(Path(a.drug_sets).read_text()).get("drug_sets", {}) if Path(a.drug_sets).exists() else {}
    wanted = set(drug_sets.get(a.drug_set, [])) if a.drug_set else None

    rows, coverage = [], []
    for npi_id, grp in npi.groupby("npi_id"):
        for (cl, mid), g in grp.dropna(subset=["model_id"]).groupby(["cell_line_id", "model_id"]):
            pick = g[g["assay"] == a.assay]
            pick = pick if len(pick) else g
            ni = float(pick["inhibition"].iloc[0])
            assay = pick["assay"].iloc[0]
            found_any = False
            for sname, auc in sources.items():
                if mid not in auc.index:
                    continue
                s = auc.loc[mid].dropna()
                if wanted:
                    s = s[[d for d in s.index if d in wanted]]
                for drug, v in s.items():
                    di = float(auc_to_inhibition(v))
                    rows.append({"npi_id": npi_id, "cell_line_id": cl,
                                 "cell_line": name_of.get(mid, cl), "drug": drug,
                                 "npi_assay": assay, "npi_inhibition": ni, "drug_source": sname,
                                 "drug_auc": float(v), "drug_inhibition": di,
                                 "expected_combined_inhibition": float(bliss_expected(ni, di))})
                    found_any = True
            coverage.append({"npi_id": npi_id, "cell_line_id": cl, "npi_inhibition": ni,
                             "npi_assay": assay, "drug_response_available": found_any})
    t = pd.DataFrame(rows)
    cov = pd.DataFrame(coverage)
    if t.empty:
        print("no (NPI, cell line) pair has both a curated NPI monotherapy value and drug response data",
              file=sys.stderr)
        cov.to_csv(out / "coverage.tsv", sep="\t", index=False)
        return 1
    t = t.sort_values("expected_combined_inhibition", ascending=False).reset_index(drop=True)
    t.to_csv(out / "ranked.tsv", sep="\t", index=False)
    cov.to_csv(out / "coverage.tsv", sep="\t", index=False)

    summary = {
        "quantity": "expected combined inhibition (total kill), Bliss independence of measured single-agent effects",
        "not_a_synergy_prediction": "measured efficacy and Bliss synergy correlate 0.060 in DrugComb",
        "validation": "within-cell-line Spearman 0.752 for ranking measured combination efficacy on "
                      "739,964 DrugComb drug-pair experiments (docs/efficacy_from_monotherapy.md)",
        "n_rows": len(t), "n_npis": int(t["npi_id"].nunique()),
        "n_cell_lines_ranked": int(t["cell_line_id"].nunique()),
        "n_cell_lines_with_npi_monotherapy": int(cov["cell_line_id"].nunique()),
        "cell_lines_ranked": sorted(t["cell_line"].unique().tolist()),
        "npi_assay_used": a.assay,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    print(f"NPI monotherapy curated for {summary['n_cell_lines_with_npi_monotherapy']} cell line(s); "
          f"ranked {len(t)} NPI x drug rows in {summary['n_cell_lines_ranked']} line(s)")
    show = ["npi_id", "cell_line", "drug", "npi_inhibition", "drug_inhibition",
            "expected_combined_inhibition", "drug_source"]
    print(t[show].head(20).to_string(index=False))
    print("\nper cell line, best-scoring drug:")
    print(t.loc[t.groupby(["npi_id", "cell_line"])["expected_combined_inhibition"].idxmax(), show].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
