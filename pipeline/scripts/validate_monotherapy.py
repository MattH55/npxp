"""What the reversal score DOES predict: single-agent sensitivity per cell line.

In the cell-matched synergy test, `mono_reversal` correlated with synergy but with
the sign flipping between metrics (negative for Bliss and ZIP, positive for HSA).
That is the fingerprint of a monotherapy-potency confound rather than synergy:
when both drugs are potent in a line, Bliss excess falls toward its ceiling while
HSA excess rises. So test the potency claim head-on.

For each (drug, cell line): reversal = -cos(s_P, s_drug) on the line's focus
genes, where s_P is the line's lineage-centred expression z-score and s_drug its
LINCS signature. Prediction: higher reversal means lower AUC (more sensitive).

The test is the residual Spearman after removing each drug's and each cell
line's mean AUC, so neither "this drug kills everything" nor "this line is
fragile" can produce the correlation.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from npi_pharma.cancer.programs import center_within, load_expression, load_gdsc_auc, load_prism_auc
from npi_pharma.store import load_signatures


def residual_spearman(t: pd.DataFrame, x: str, y: str, g1: str, g2: str) -> tuple[float, float, int]:
    v = t[y]
    resid = v - t.groupby(g1)[y].transform("mean") - t.groupby(g2)[y].transform("mean") + v.mean()
    return (float(t[x].corr(v, method="spearman")), float(t[x].corr(resid, method="spearman")), len(t))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--lincs", required=True, help="consensus drug signature parquet")
    ap.add_argument("--expression", default="data/raw/depmap/OmicsExpressionProteinCodingGenesTPMLogp1.csv")
    ap.add_argument("--models", default="data/raw/depmap/Model.csv")
    ap.add_argument("--prism", default="data/raw/depmap/prism19q4_secondary_dose_response.csv")
    ap.add_argument("--gdsc", nargs="*", default=["data/raw/depmap/GDSC1_fitted_dose_response_25Feb20.csv",
                                                  "data/raw/depmap/GDSC2_fitted_dose_response_25Feb20.csv"])
    ap.add_argument("--focus-top-k", type=int, default=500)
    ap.add_argument("--permutations", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="out/drugcomb_validation")
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    sigs = pd.DataFrame({s.sig_id: s.as_series() for s in load_signatures(a.lincs)})
    models = pd.read_csv(a.models)
    lineage = models.set_index("ModelID")["OncotreeLineage"]
    expr = load_expression(a.expression, list(sigs.index))
    genes = sigs.index.intersection(expr.columns)
    sigs = sigs.loc[genes]
    xc = center_within(expr[genes], lineage)
    P = ((xc - xc.mean()) / xc.std(ddof=0)).T                    # genes x cell lines
    unit = sigs / np.sqrt((sigs ** 2).sum())
    k = min(a.focus_top_k, len(genes))
    print(f"{len(genes)} genes; {sigs.shape[1]} drugs; {P.shape[1]} cell lines", file=sys.stderr)

    # reversal for every (cell line, drug)
    pm = P.to_numpy()
    foc = np.argsort(-np.abs(pm), axis=0)[:k]                    # k x lines
    pf = np.take_along_axis(pm, foc, axis=0)
    pn = np.linalg.norm(pf, axis=0)
    rv = {}
    for drug in unit.columns:
        v = unit[drug].to_numpy()[foc]
        rv[drug] = -(pf * v).sum(0) / (pn * np.linalg.norm(v, axis=0))
    R = pd.DataFrame(rv, index=P.columns)                        # cell lines x drugs
    R.to_pickle(out / "reversal_cellline_by_drug.pkl")

    sources = {"PRISM": load_prism_auc(a.prism)}
    for p in a.gdsc:
        sources[Path(p).name.split("_")[0]] = load_gdsc_auc(p, models)

    res = {"n_genes": int(len(genes)), "n_drugs": int(R.shape[1]), "n_cell_lines": int(R.shape[0]),
           "focus_top_k": a.focus_top_k}
    rng = np.random.default_rng(a.seed)
    for name, auc in sources.items():
        drugs = [d for d in R.columns if d in auc.columns]
        lines = R.index.intersection(auc.index)
        if not drugs or not len(lines):
            continue
        t = (auc.loc[lines, drugs].stack().rename("auc").reset_index()
             .rename(columns={"level_0": "cell", "level_1": "drug"}))
        t.columns = ["cell", "drug", "auc"]
        t["reversal"] = [R.at[c, d] for c, d in zip(t["cell"], t["drug"])]
        raw, resid, n = residual_spearman(t, "reversal", "auc", "drug", "cell")
        v = t["auc"]
        r0 = (v - t.groupby("drug")["auc"].transform("mean") - t.groupby("cell")["auc"].transform("mean") + v.mean())
        null = np.array([float(pd.Series(rng.permutation(t["reversal"].to_numpy())).corr(r0, method="spearman"))
                         for _ in range(a.permutations)])
        res[name] = {"n_observations": n, "n_drugs": len(drugs), "n_cell_lines": int(len(lines)),
                     "spearman_raw": raw, "spearman_residual": resid,
                     "residual_permutation_p": float((1 + (np.abs(null) >= abs(resid)).sum()) / (1 + a.permutations)),
                     "residual_null_sd": float(null.std())}
        # per-drug residual correlation, to see whether it is a few drugs or general
        per = []
        for drug, g in t.groupby("drug"):
            if len(g) < 100:
                continue
            gv = g["auc"] - g["auc"].mean() - t.set_index("cell").loc[g["cell"], "auc"].groupby(level=0).mean().reindex(g["cell"]).to_numpy() + v.mean()
            per.append(g["reversal"].corr(pd.Series(gv.to_numpy(), index=g.index), method="spearman"))
        per = pd.Series(per).dropna()
        res[name]["per_drug_residual_spearman"] = {"n_drugs": int(len(per)), "median": float(per.median()),
                                                  "frac_negative": float((per < 0).mean())}
        t.to_csv(out / f"monotherapy_{name}.tsv", sep="\t", index=False)
    (out / "monotherapy_summary.json").write_text(json.dumps(res, indent=2, default=float))
    for name in sources:
        if name in res:
            e = res[name]
            print(f"{name:6s} n={e['n_observations']:7d} ({e['n_drugs']} drugs x {e['n_cell_lines']} lines)  "
                  f"rho_raw {e['spearman_raw']:+.3f}  rho_resid {e['spearman_residual']:+.3f} "
                  f"(p={e['residual_permutation_p']:.4f}, null sd {e['residual_null_sd']:.3f})  "
                  f"per-drug median {e['per_drug_residual_spearman']['median']:+.3f}, "
                  f"{e['per_drug_residual_spearman']['frac_negative']:.0%} negative")
    print("negative rho = higher reversal -> lower AUC = more sensitive (the predicted direction)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
