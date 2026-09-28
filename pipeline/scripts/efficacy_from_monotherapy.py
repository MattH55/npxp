"""Apply the literature's finding: rank combinations by predicted EFFICACY, not synergy.

The DREAM challenge and the review literature agree that each agent's own response
is the strongest single predictor of combination outcome, and that predicting
excess-over-additivity for unseen drugs remains unsolved. This project's own
tests agree: signature composition reaches |rho| <= 0.02 on measured synergy
(docs/validation_drugcomb.md).

In DrugComb, combination efficacy (CSS, how much the combination actually
inhibits) correlates 0.69 with the agents' own relative inhibition and only 0.06
with Bliss synergy. Those are nearly independent quantities, and it is efficacy
that decides whether a combination kills the cells.

So this script asks the question that can be answered: given both agents'
monotherapy responses in a cell line, how well can the combination's efficacy in
that cell line be predicted? Validation is leave-one-cell-line-out (GroupKFold by
cell line), i.e. every test cell line is unseen, which is the regime an
application needs.

Three predictors are compared:
  pair_mean       the pair's mean efficacy in the OTHER cell lines. No knowledge
                  of this cell line at all; the "just use the average" baseline.
  monotherapy     both agents' relative inhibition in this cell line, and simple
                  combinations of the two. No pair identity.
  both            monotherapy features plus the pair mean.

Also reported, for contrast, is the same comparison with Bliss synergy as the
target, to show the gap between the answerable and unanswerable question.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import GroupKFold

FEATURES = ["ri_row", "ri_col", "ri_min", "ri_max", "ri_mean", "ri_diff"]


def add_features(d: pd.DataFrame) -> pd.DataFrame:
    d = d.copy()
    d["ri_min"] = d[["ri_row", "ri_col"]].min(axis=1)
    d["ri_max"] = d[["ri_row", "ri_col"]].max(axis=1)
    d["ri_mean"] = d[["ri_row", "ri_col"]].mean(axis=1)
    d["ri_diff"] = (d["ri_row"] - d["ri_col"]).abs()
    return d


def evaluate(d: pd.DataFrame, target: str, n_splits: int, seed: int) -> dict:
    """Leave-cell-lines-out scores for the three predictors."""
    groups = d["cell_line_name"].to_numpy()
    y = d[target].to_numpy()
    preds = {k: np.full(len(d), np.nan) for k in ("pair_mean", "monotherapy", "both")}
    for tr, te in GroupKFold(n_splits=n_splits).split(d, y, groups):
        train, test = d.iloc[tr], d.iloc[te]
        # pair_mean: the pair's mean target over training cell lines (global mean if unseen pair)
        pm_train = train.groupby("pair")[target].mean()
        gmean = float(train[target].mean())
        pm_te = test["pair"].map(pm_train).fillna(gmean).to_numpy()
        pm_tr = train["pair"].map(pm_train).fillna(gmean).to_numpy()
        preds["pair_mean"][te] = pm_te
        for name, cols, extra_tr, extra_te in (
            ("monotherapy", FEATURES, None, None),
            ("both", FEATURES, pm_tr, pm_te),
        ):
            xtr, xte = train[cols].to_numpy(), test[cols].to_numpy()
            if extra_tr is not None:
                xtr = np.column_stack([xtr, extra_tr])
                xte = np.column_stack([xte, extra_te])
            m = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.1, random_state=seed)
            m.fit(xtr, y[tr])
            preds[name][te] = m.predict(xte)
    out = {}
    for k, p in preds.items():
        ok = ~np.isnan(p)
        out[k] = {
            "spearman": float(pd.Series(p[ok]).corr(pd.Series(y[ok]), method="spearman")),
            "pearson": float(pd.Series(p[ok]).corr(pd.Series(y[ok]))),
            "rmse": float(np.sqrt(np.mean((p[ok] - y[ok]) ** 2))),
        }
        # per-cell-line ranking quality: can it order combinations within a new line?
        per = [pd.Series(p[g.index]).corr(pd.Series(y[g.index]), method="spearman")
               for _, g in d.reset_index(drop=True).assign(_p=p).groupby("cell_line_name") if len(g) >= 20]
        per = pd.Series(per).dropna()
        out[k]["within_cell_line_spearman_median"] = float(per.median())
        out[k]["n_cell_lines_scored"] = int(len(per))
    out["_target_sd"] = float(np.std(y))
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--drugcomb", required=True)
    ap.add_argument("--splits", type=int, default=5)
    ap.add_argument("--min-rows-per-cell-line", type=int, default=50)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="out/efficacy_from_monotherapy")
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    cols = ["drug_row", "drug_col", "cell_line_name", "study_name",
            "ri_row", "ri_col", "css_row", "css_col", "synergy_bliss"]
    d = pd.read_csv(a.drugcomb, usecols=cols, low_memory=False)
    d = d[d["drug_col"].notna() & (d["drug_col"] != "NULL")]
    for c in cols[4:]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.dropna(subset=["ri_row", "ri_col", "css_row", "css_col", "synergy_bliss"])
    d["efficacy_css"] = d[["css_row", "css_col"]].mean(axis=1)
    d["pair"] = [" + ".join(sorted(p)) for p in zip(d["drug_row"], d["drug_col"])]
    keep = d["cell_line_name"].map(d["cell_line_name"].value_counts()) >= a.min_rows_per_cell_line
    d = add_features(d[keep]).reset_index(drop=True)
    print(f"{len(d)} rows, {d['pair'].nunique()} pairs, {d['cell_line_name'].nunique()} cell lines",
          file=sys.stderr)

    res = {"n_rows": len(d), "n_pairs": int(d["pair"].nunique()),
           "n_cell_lines": int(d["cell_line_name"].nunique()), "splits": a.splits,
           "validation": "GroupKFold by cell line (every test cell line unseen)",
           "correlation_efficacy_vs_synergy": float(d["efficacy_css"].corr(d["synergy_bliss"])),
           "correlation_efficacy_vs_monotherapy": float(d["efficacy_css"].corr(d["ri_mean"]))}
    for target in ("efficacy_css", "synergy_bliss"):
        res[target] = evaluate(d, target, a.splits, a.seed)
    (out / "summary.json").write_text(json.dumps(res, indent=2, default=float))

    print(f"rows {res['n_rows']}, pairs {res['n_pairs']}, cell lines {res['n_cell_lines']}; "
          f"corr(efficacy, synergy) = {res['correlation_efficacy_vs_synergy']:.3f}")
    for target in ("efficacy_css", "synergy_bliss"):
        print(f"\n{target} (sd {res[target]['_target_sd']:.2f}), leave-cell-lines-out:")
        for k in ("pair_mean", "monotherapy", "both"):
            e = res[target][k]
            print(f"  {k:12s} spearman {e['spearman']:+.3f}  pearson {e['pearson']:+.3f}  "
                  f"rmse {e['rmse']:.2f}  within-cell-line median rho {e['within_cell_line_spearman_median']:+.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
