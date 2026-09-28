"""Audit of the Phase 7 modifier x drug x cell-line predictions.

Phase 7 (pathway-feature GBT, Cancer Resource/synlethality) passed its RMSE gate
against Bliss independence. This script asks whether its predictions can rank
combinations, using checks it was not trained or tuned on:

1. Negative controls. Phase 3 registered three modifiers with the weakest
   mimetic anchoring as negative controls; they should not outscore the rest.
2. Error bars. How many predictions exceed an approximate 90% interval
   half-width, taken as 1.645 x the drug-out RMSE, because phase7 does not store
   the conformal quantile.
3. Variance decomposition. Main effects of modifier, drug and cell line vs. the
   interaction residual.
4. Same-trial arms. GSE85620 placebo and cold-water-immersion arms share the
   same 10-week strength programme.
5. Directional anchor check. The three curated tier-1 heat x cisplatin enhancement
   ratios (tests/test_seed.py REAL_COMBINED_EFFECT_METRICS) that fall inside
   the Phase 7 grid. The metric differs (enhancement ratio vs Bliss excess), so
   only direction and rank are compared: heat should score positive and above
   other modifiers for cisplatin in those cell lines.

Input is a checkout of MattH55/Physiological-Fitness-Landscape (the
`Cancer Resource` folder, commit 7f72fbb or later):

    python analysis/phase7_audit/audit_phase7.py --cancer-resource "<path>/Cancer Resource"
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

Y = "predicted_bliss_excess"
CONTROLS = ["gse153830_t47d_bhb25mm", "gse28016_muscle_fast40h_vs_fed", "gse85620_cwi_post_vs_pre"]
HEAT_IN_VITRO = [
    "gse48398_mcf7_heat45c30min", "gse48398_mda231_heat45c30min", "gse48398_mda468_heat45c30min",
    "gse48398_mcf10a_heat45c30min", "gse10043_u937_mildhyperthermia41c30min", "gse75127_hsc3_hyperthermia44c90min",
]
HEAT_IN_VIVO = ["gse82323_muscle_heat", "gse12474_muscle_heat_sheet10w", "gse90763_pbmc_sauna_15min_after"]
# tests/test_seed.py REAL_COMBINED_EFFECT_METRICS entries whose drug and cell line are in the Phase 7 grid.
ANCHORS = [
    ("MOD-HT-42.8C-30M-SIMUL", "cisplatin", "HELA_CERVIX", 3.10, "Kusumoto 1993, PMID 8347479"),
    ("MOD-HT-43C-60M-HIPEC", "cisplatin", "RKO_LARGE_INTESTINE", 3.5, "Helderman 2020"),
    ("MOD-HT-43C-60M-HIPEC", "cisplatin", "HCT116_LARGE_INTESTINE", 2.8, "Helderman 2020"),
]


def auc(pos: np.ndarray, neg: np.ndarray) -> float:
    r = pd.Series(np.r_[pos, neg]).rank().to_numpy()
    return float((r[: len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cancer-resource", required=True, help="path to the 'Cancer Resource' folder")
    ap.add_argument("--out", help="write the results as JSON here")
    a = ap.parse_args(argv)
    root = Path(a.cancer_resource)
    d = json.loads((root / "data/lincs/phase7_predictions.json").read_text())
    p = pd.DataFrame(d["predictions"])
    res: dict = {"n_predictions": len(p), "drug_out_rmse": d["gbt_rmse_drug_out"],
                 "coverage_90_drug_out": d["coverage_90_drug_out"]}

    # 1. negative controls
    is_ctrl = p["modifier_signature_id"].isin(CONTROLS)
    by_mod = p.groupby("modifier_signature_id")[Y].mean().sort_values(ascending=False)
    rank = {m: int(i) + 1 for i, m in enumerate(by_mod.index)}
    res["controls"] = {
        "auc_real_modifier_above_control": auc(p.loc[~is_ctrl, Y].to_numpy(), p.loc[is_ctrl, Y].to_numpy()),
        "control_ranks_of_33": {c: rank[c] for c in CONTROLS},
    }

    # 2. error bars
    half = 1.645 * d["gbt_rmse_drug_out"]
    res["approx_90_half_width"] = half
    res["n_above_half_width"] = int((p[Y] > half).sum())

    # 3. variance decomposition (balanced grid)
    mu, tot = p[Y].mean(), ((p[Y] - p[Y].mean()) ** 2).sum()
    shares, resid = {}, p[Y] - mu
    for f in ("modifier_signature_id", "drug_id", "cell_line_id"):
        m = p.groupby(f)[Y].transform("mean")
        shares[f] = float(((m - mu) ** 2).sum() / tot)
        resid = resid - (m - mu)
    shares["interaction_residual"] = float((resid ** 2).sum() / tot)
    res["variance_share"] = shares

    # 4. same-trial arms
    res["gse85620_arms"] = {m: {"mean": float(by_mod[m]), "rank_of_33": rank[m]}
                            for m in ("gse85620_placebo_post_vs_pre", "gse85620_cwi_post_vs_pre")}

    # 5. directional anchor check
    rows = []
    for mod, drug, cl, er, src in ANCHORS:
        s = p[(p.drug_id == drug) & (p.cell_line_id == cl)].set_index("modifier_signature_id")[Y]
        r = s.rank(ascending=False)
        for grp, ms in (("heat_in_vitro", HEAT_IN_VITRO), ("heat_in_vivo", HEAT_IN_VIVO)):
            rows.append({"anchor": mod, "drug": drug, "cell_line": cl, "enhancement_ratio": er, "source": src,
                         "group": grp, "mean_pred": float(s[ms].mean()), "n_positive": int((s[ms] > 0).sum()),
                         "n": len(ms), "mean_rank_of_33": float(r[ms].mean()),
                         "all_modifiers_mean": float(s.mean())})
    res["anchor_check"] = rows

    top = (p.groupby(["modifier_signature_id", "drug_id"])[Y].agg(["mean", "min", "max"])
           .sort_values("mean", ascending=False).head(10).round(3))
    print(f"predictions: {len(p)}; approx 90% half-width {half:.3f}; above it: {res['n_above_half_width']}")
    print(f"controls: AUC real > control {res['controls']['auc_real_modifier_above_control']:.3f}; "
          f"ranks {res['controls']['control_ranks_of_33']}")
    print("variance share: " + ", ".join(f"{k} {v:.2f}" for k, v in shares.items()))
    print(f"GSE85620 arms: {res['gse85620_arms']}")
    print(pd.DataFrame(rows)[["cell_line", "enhancement_ratio", "group", "mean_pred", "n_positive", "n",
                              "mean_rank_of_33", "all_modifiers_mean"]].round(3).to_string(index=False))
    print("top modifier x drug (mean over 13 cell lines):\n" + top.to_string())
    if a.out:
        res["top_modifier_drug"] = top.reset_index().to_dict(orient="records")
        Path(a.out).write_text(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
