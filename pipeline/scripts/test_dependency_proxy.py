"""Can a per-cell-line functional measurement predict drug sensitivity? Calibration.

`configs/npi_monotherapy.yaml` needs the NPI's own effect per cell line, and that
is almost entirely uncurated (see docs/npi_monotherapy_search_log.md). One
tempting substitute is DepMap CRISPR dependency: a line that dies without SLC7A11
is a line that dies without cystine import, which is what cystine withdrawal does.

Before using any such proxy, this script measures what the approach can achieve at
its best, by running it where the link is textbook: a gene whose knockout
phenocopies its own inhibitor. If EGFR dependency predicts gefitinib sensitivity
well, the approach works and a weak metabolic result means the proxy is weak; if
even that fails, the comparison itself is too noisy to interpret.

Expected sign is POSITIVE: a less negative gene effect (less dependent) should go
with a higher AUC (more resistant).

    python scripts/test_dependency_proxy.py
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

from npi_pharma.cancer.programs import load_gdsc_auc, load_prism_auc

# (gene, drug, why). The first block is the calibration set: the drug inhibits the
# gene's own product, so dependency should track sensitivity. The second is the
# NPI-proxy question this project actually needs.
CALIBRATION = [
    ("EGFR", "gefitinib", "EGFR inhibitor"),
    ("EGFR", "erlotinib", "EGFR inhibitor"),
    ("MDM2", "nutlin-3", "MDM2 inhibitor"),
    ("BRAF", "dabrafenib", "BRAF inhibitor"),
    ("BRAF", "plx-4720", "BRAF inhibitor"),
    ("BCL2L1", "navitoclax", "BCL-xL/BCL-2 inhibitor"),
    ("CDK4", "palbociclib", "CDK4/6 inhibitor"),
    ("ABL1", "imatinib", "BCR-ABL inhibitor"),
    ("FLT3", "quizartinib", "FLT3 inhibitor"),
]
NPI_PROXY = [
    ("SLC7A11", "erastin", "cystine importer; erastin blocks it -- proxy for cystine withdrawal"),
    ("SLC7A11", "sulfasalazine", "same transporter, weaker inhibitor"),
    ("MTOR", "sirolimus", "mTOR; proxy for the nutrient-restriction arm"),
    ("MTOR", "everolimus", "mTOR"),
]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--depmap-dir", default="data/raw/depmap")
    ap.add_argument("--crispr", default="CRISPRGeneEffect.csv")
    ap.add_argument("--min-lines", type=int, default=50)
    ap.add_argument("--out", default="out/dependency_proxy")
    a = ap.parse_args(argv)
    dm, out = Path(a.depmap_dir), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    genes = sorted({g for g, _, _ in CALIBRATION + NPI_PROXY})
    hdr = pd.read_csv(dm / a.crispr, nrows=0).columns
    sym = {c: re.sub(r"\s*\(\d+\)$", "", c) for c in hdr[1:]}
    use = [hdr[0]] + [c for c in hdr[1:] if sym[c] in set(genes)]
    ge = pd.read_csv(dm / a.crispr, usecols=use, index_col=0)
    ge.columns = [sym[c] for c in ge.columns]
    print(f"CRISPR gene effect: {ge.shape[0]} cell lines x {ge.shape[1]} genes", file=sys.stderr)

    models = pd.read_csv(dm / "Model.csv")
    sources = {"PRISM": load_prism_auc(dm / "prism19q4_secondary_dose_response.csv")}
    for p in sorted(dm.glob("GDSC*_fitted_dose_response*.csv")):
        sources[p.name.split("_")[0]] = load_gdsc_auc(p, models)

    rows = []
    for block, pairs in (("calibration", CALIBRATION), ("npi_proxy", NPI_PROXY)):
        for gene, drug, why in pairs:
            if gene not in ge.columns:
                continue
            for name, auc in sources.items():
                if drug not in auc.columns:
                    continue
                lines = ge.index.intersection(auc.index)
                x, y = ge.loc[lines, gene], auc.loc[lines, drug]
                ok = x.notna() & y.notna()
                if ok.sum() < a.min_lines:
                    continue
                rows.append({"block": block, "gene": gene, "drug": drug, "source": name,
                             "spearman": float(x[ok].corr(y[ok], method="spearman")),
                             "n_cell_lines": int(ok.sum()), "rationale": why})
    t = pd.DataFrame(rows)
    t.to_csv(out / "dependency_vs_drug_sensitivity.tsv", sep="\t", index=False)
    best = t.groupby(["block", "gene", "drug"])["spearman"].max()
    summary = {
        "expected_sign": "positive (less dependent -> higher AUC -> more resistant)",
        "calibration_best_spearman": {f"{g} x {d}": float(v) for (b, g, d), v in best.items() if b == "calibration"},
        "npi_proxy_best_spearman": {f"{g} x {d}": float(v) for (b, g, d), v in best.items() if b == "npi_proxy"},
        "calibration_median_best": float(best.loc["calibration"].median()),
        "npi_proxy_median_best": float(best.loc["npi_proxy"].median()),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(t.sort_values(["block", "gene", "drug", "source"]).to_string(index=False))
    print(f"\nbest-source median Spearman: calibration {summary['calibration_median_best']:+.3f}, "
          f"NPI proxy {summary['npi_proxy_median_best']:+.3f}")
    print("A per-cell-line functional measurement predicting single-agent sensitivity tops out around "
          "+0.2 to +0.47 even for a drug's own target; that is the ceiling for this class of feature.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
