"""Does the npi_pharma pair score predict which CELL LINES show synergy?

The pair-level test (scripts/validate_reversal_drugcomb.py) asks whether a pair's
score predicts its median synergy. This asks the question the tool is actually
used for: given a cell line, does the score pick the combinations that work in
THAT line?

For each (drug A, drug B, cell line C) in DrugComb with DepMap expression for C
and LINCS signatures for A and B, compute the npi_pharma scores with the cell
line in the patient slot: s_P is C's lineage-centred expression z-score, and
reversal, complementarity and mono_reversal follow npi_pharma.interact.score.

The test is a within-pair, within-cell-line residual: synergy minus the pair's
mean minus the cell line's mean. A pair that synergises everywhere, or a cell
line that shows high synergy with everything, must not count as a hit.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from npi_pharma.cancer.programs import center_within
from npi_pharma.store import load_signatures

SYNERGY = ["synergy_bliss", "synergy_loewe", "synergy_zip", "synergy_hsa"]
NORM = lambda s: re.sub(r"[^a-z0-9]", "", str(s).lower())


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--drugcomb", required=True)
    ap.add_argument("--lincs", required=True)
    ap.add_argument("--expression", default="data/raw/depmap/OmicsExpressionProteinCodingGenesTPMLogp1.csv")
    ap.add_argument("--models", default="data/raw/depmap/Model.csv")
    ap.add_argument("--focus-top-k", type=int, default=500)
    ap.add_argument("--permutations", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="out/drugcomb_validation")
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    sigs = pd.DataFrame({s.sig_id: s.as_series() for s in load_signatures(a.lincs)})
    models = pd.read_csv(a.models)
    lineage = models.set_index("ModelID")["OncotreeLineage"]
    from npi_pharma.cancer.programs import load_expression

    expr = load_expression(a.expression, list(sigs.index))
    genes = sigs.index.intersection(expr.columns)
    sigs = sigs.loc[genes]
    # s_P: each line's expression z against other lines of its lineage
    xc = center_within(expr[genes], lineage)
    P = ((xc - xc.mean()) / xc.std(ddof=0)).T  # genes x cell lines
    print(f"{len(genes)} genes; {sigs.shape[1]} signatures; {P.shape[1]} cell lines", file=sys.stderr)

    d = pd.read_csv(a.drugcomb, usecols=["drug_row", "drug_col", "cell_line_name", "study_name", *SYNERGY],
                    low_memory=False)
    d = d[d["drug_col"].notna() & (d["drug_col"] != "NULL")]
    for c in SYNERGY:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.dropna(subset=["synergy_bliss"])
    nl = {NORM(c): c for c in sigs.columns}
    d["A"] = d["drug_row"].map(lambda x: nl.get(NORM(x)))
    d["B"] = d["drug_col"].map(lambda x: nl.get(NORM(x)))
    mm = {NORM(r.StrippedCellLineName): r.ModelID for r in models.itertuples()
          if isinstance(r.StrippedCellLineName, str)}
    d["ModelID"] = d["cell_line_name"].map(lambda c: mm.get(NORM(c)))
    d = d.dropna(subset=["A", "B", "ModelID"])
    d = d[d["ModelID"].isin(P.columns) & (d["A"] != d["B"])]
    d["pair"] = [" + ".join(sorted(p)) for p in zip(d["A"], d["B"])]
    d = d.groupby(["pair", "A", "B", "ModelID"], as_index=False)[SYNERGY].median()
    print(f"{len(d)} (pair, cell line) observations; {d['pair'].nunique()} pairs, "
          f"{d['ModelID'].nunique()} cell lines", file=sys.stderr)

    # npi_pharma scores, per (pair, cell line)
    unit = sigs / np.sqrt((sigs ** 2).sum())
    rows = []
    for (pair, A, B), grp in d.groupby(["pair", "A", "B"]):
        sa, sb = unit[A].to_numpy(), unit[B].to_numpy()
        combo = sa + sb
        p = P[grp["ModelID"].to_numpy()].to_numpy()          # genes x lines
        k = min(a.focus_top_k, len(genes))
        foc = np.argsort(-np.abs(p), axis=0)[:k]             # per-line focus genes
        pf = np.take_along_axis(p, foc, axis=0)          # k x lines
        pn = np.linalg.norm(pf, axis=0)

        def rv(v: np.ndarray) -> np.ndarray:
            """-cos(s_P, v) on each line's own focus genes; v is 1-D over genes."""
            vf = v[foc]                                   # k x lines
            return -(pf * vf).sum(0) / (pn * np.linalg.norm(vf, axis=0))

        r_a, r_b, r_c = rv(sa), rv(sb), rv(combo)
        best = np.maximum(r_a, r_b)
        rows.append(pd.DataFrame({
            "pair": pair, "ModelID": grp["ModelID"].to_numpy(),
            "mono_reversal": (r_a + r_b) / 2,
            "complementarity_gain": r_c - best,
            "complementarity": (r_c - best) / np.maximum(1 - best, 1e-6),
            "reverse_combo": r_c,
            **{c: grp[c].to_numpy() for c in SYNERGY},
        }))
    t = pd.concat(rows, ignore_index=True)
    t.to_csv(out / "cellwise_scores.tsv", sep="\t", index=False)

    feats = ["complementarity", "complementarity_gain", "mono_reversal", "reverse_combo"]
    res = {"n_observations": len(t), "n_pairs": int(t["pair"].nunique()),
           "n_cell_lines": int(t["ModelID"].nunique()), "n_genes": int(len(genes)),
           "focus_top_k": a.focus_top_k}
    rng = np.random.default_rng(a.seed)
    for metric in SYNERGY:
        y = t[metric]
        resid = (y - t.groupby("pair")[metric].transform("mean")
                 - t.groupby("ModelID")[metric].transform("mean") + y.mean()).to_numpy()
        e = {}
        for f in feats:
            x = t[f].to_numpy()
            e[f] = {"spearman_raw": float(pd.Series(x).corr(y, method="spearman")),
                    "spearman_residual": float(pd.Series(x).corr(pd.Series(resid), method="spearman"))}
        # permutation null for the primary feature, shuffling cell lines within each pair
        prim = t["complementarity"].to_numpy()
        obs = e["complementarity"]["spearman_residual"]
        null = []
        gi = t.groupby("pair").indices
        for _ in range(a.permutations):
            sh = prim.copy()
            for idx in gi.values():
                sh[idx] = rng.permutation(sh[idx])
            null.append(float(pd.Series(sh).corr(pd.Series(resid), method="spearman")))
        null = np.array(null)
        e["complementarity"]["residual_permutation_p"] = float((1 + (np.abs(null) >= abs(obs)).sum()) / (1 + a.permutations))
        e["complementarity"]["residual_null_sd"] = float(null.std())
        res[metric] = e
    (out / "cellwise_summary.json").write_text(json.dumps(res, indent=2, default=float))
    print(f"observations: {res['n_observations']}  pairs: {res['n_pairs']}  cell lines: {res['n_cell_lines']}")
    for metric in SYNERGY:
        e = res[metric]
        print(f"{metric:15s} " + "  ".join(
            f"{f}: raw {e[f]['spearman_raw']:+.3f} resid {e[f]['spearman_residual']:+.3f}" for f in feats)
            + f"  (complementarity p={e['complementarity']['residual_permutation_p']:.4f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
