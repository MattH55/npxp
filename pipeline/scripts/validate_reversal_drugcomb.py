"""Does resistance-program reversal predict measured synergy? Test on DrugComb drug x drug.

The NPI x drug predictor (npi_pharma.cancer) assumes: agent X sensitises cells to
drug Y when X's expression signature moves cells toward the baseline state of
lines that are sensitive to Y. For NPI x drug there are almost no measured
outcomes to check that against. For drug x drug there are 740k, in DrugComb.

So score drug pairs by exactly the same rule, using LINCS L1000 signatures for
one side and the PRISM/DepMap resistance program for the other:

    score(A, B) = mean( -corr(s_A, r_B), -corr(s_B, r_A) )

and test it against measured synergy. Two tests:

  raw       Spearman of the score with pair-level median synergy.
  residual  Same, after removing each drug's own mean synergy (its main effect),
            i.e. does the score explain what is specific to the pair? A broadly
            cytotoxic drug synergises with many partners and also has a strong
            signature, which the raw test cannot separate.

Baseline for comparison: |cos| between the two LINCS signatures (the
monotherapy-correlation feature, eLife 2020;9:e52707).

Inputs: out/cancer_npi_drug/programs_cache.pkl (from scripts/cancer_npi_drug.py),
a LINCS signature parquet for the drugs, and DrugComb summary_v_1_5.csv.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from npi_pharma.store import load_signatures

SYNERGY = ["synergy_bliss", "synergy_loewe", "synergy_zip", "synergy_hsa"]
NORM = lambda s: re.sub(r"[^a-z0-9]", "", str(s).lower())


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return float(pd.Series(x).corr(pd.Series(y), method="spearman"))


def auc(score: np.ndarray, label: np.ndarray) -> float:
    r = pd.Series(score).rank().to_numpy()
    n1, n0 = int(label.sum()), int((~label.astype(bool)).sum())
    if not n1 or not n0:
        return float("nan")
    return float((r[label.astype(bool)].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--drugcomb", required=True, help="DrugComb summary_v_1_5.csv")
    ap.add_argument("--programs", default="out/cancer_npi_drug/programs_cache.pkl")
    ap.add_argument("--lincs", required=True, help="drug signature parquet from ingest-lincs")
    ap.add_argument("--min-cell-lines", type=int, default=3, help="cell lines a pair needs for a median")
    ap.add_argument("--synergy-threshold", type=float, default=10.0, help="|synergy| for the AUC classes")
    ap.add_argument("--permutations", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="out/drugcomb_validation")
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    programs, specific, _sources, _n, _rel = pd.read_pickle(a.programs)
    sigs = {s.sig_id: s.as_series() for s in load_signatures(a.lincs)}
    S = pd.DataFrame(sigs)
    genes = S.index.intersection(specific.index)
    Sz = S.loc[genes].apply(lambda c: (c - c.mean()) / c.std(ddof=0))
    Rz = specific.loc[genes].apply(lambda c: (c - c.mean()) / c.std(ddof=0))
    rev = -(Sz.T.fillna(0) @ Rz.fillna(0)) / len(genes)          # signature x program
    cos = (Sz.T.fillna(0) @ Sz.fillna(0)) / len(genes)           # signature x signature
    print(f"{len(genes)} shared genes; {Sz.shape[1]} signatures; {Rz.shape[1]} programs", file=sys.stderr)

    d = pd.read_csv(a.drugcomb, usecols=["drug_row", "drug_col", "cell_line_name", "study_name", *SYNERGY],
                    low_memory=False)
    d = d[d["drug_col"].notna() & (d["drug_col"] != "NULL")]
    for c in SYNERGY:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.dropna(subset=["synergy_bliss"])
    nl = {NORM(c): c for c in Sz.columns}
    npr = {NORM(c): c for c in Rz.columns}
    for side, col in (("A", "drug_row"), ("B", "drug_col")):
        n = d[col].map(NORM)
        d[f"{side}_sig"] = n.map(nl)
        d[f"{side}_prog"] = n.map(npr)
    # A pair is scorable if at least one direction exists: a signature for one
    # drug and a program for the other. Both-direction pairs average the two.
    d["dir_ab"] = d["A_sig"].notna() & d["B_prog"].notna()
    d["dir_ba"] = d["B_sig"].notna() & d["A_prog"].notna()
    d = d[d["dir_ab"] | d["dir_ba"]]
    d["key_a"] = d["A_sig"].fillna(d["A_prog"])
    d["key_b"] = d["B_sig"].fillna(d["B_prog"])
    print(f"{len(d)} measured combinations scorable in at least one direction "
          f"({int((d.dir_ab & d.dir_ba).sum())} both)", file=sys.stderr)

    d["pair"] = [tuple(sorted(p)) for p in zip(d["key_a"], d["key_b"])]
    pairs = d.groupby("pair").agg(n_cell_lines=("cell_line_name", "nunique"),
                                  **{c: (c, "median") for c in SYNERGY}).reset_index()
    pairs = pairs[pairs["n_cell_lines"] >= a.min_cell_lines]
    pairs = pairs[[x != y for x, y in pairs["pair"]]]
    A = pairs["pair"].str[0]
    B = pairs["pair"].str[1]

    def pair_score(x: str, y: str) -> tuple[float, int]:
        vals = []
        for u, v in ((x, y), (y, x)):
            if u in rev.index and NORM(v) in npr:
                vals.append(rev.at[u, npr[NORM(v)]])
        return (float(np.mean(vals)) if vals else np.nan), len(vals)

    sc = [pair_score(x, y) for x, y in zip(A, B)]
    pairs["score"] = [v for v, _ in sc]
    pairs["n_directions"] = [n for _, n in sc]
    pairs["abs_cos_signatures"] = [abs(cos.at[x, y]) if x in cos.index and y in cos.columns else np.nan
                                   for x, y in zip(A, B)]
    pairs = pairs.dropna(subset=["score"])
    pairs["drug_a"], pairs["drug_b"] = A, B

    res: dict = {"n_pairs_both_directions": int((pairs["n_directions"] == 2).sum()) if "n_directions" in pairs else 0,
                 "n_pairs": len(pairs), "n_drugs": int(pd.unique(pairs[["drug_a", "drug_b"]].values.ravel()).size),
                 "n_combination_rows": int(len(d)), "n_genes": int(len(genes)),
                 "min_cell_lines_per_pair": a.min_cell_lines}
    rng = np.random.default_rng(a.seed)
    for metric in SYNERGY:
        y = pairs[metric].to_numpy()
        # drug main effects, from this pair table
        long = pd.concat([pairs[["drug_a", metric]].rename(columns={"drug_a": "drug"}),
                          pairs[["drug_b", metric]].rename(columns={"drug_b": "drug"})])
        mean_by_drug = long.groupby("drug")[metric].mean()
        resid = y - mean_by_drug[pairs["drug_a"]].to_numpy() - mean_by_drug[pairs["drug_b"]].to_numpy() + y.mean()
        entry = {}
        for name, s in (("reversal_score", pairs["score"].to_numpy()), ("abs_cos_signatures", pairs["abs_cos_signatures"].to_numpy())):
            hi = y >= a.synergy_threshold
            lo = y <= -a.synergy_threshold
            keep = (hi | lo) & ~np.isnan(s)
            fin = ~np.isnan(s)
            entry[name] = {
                "spearman_raw": spearman(s[fin], y[fin]),
                "spearman_residual": spearman(s[fin], resid[fin]),
                "auc_synergy_vs_antagonism": auc(s[keep], hi[keep]),
                "n_synergy": int((hi & fin).sum()), "n_antagonism": int((lo & fin).sum()),
                "n_scored": int(fin.sum()),
            }
        null = np.array([spearman(rng.permutation(pairs["score"].to_numpy()), resid) for _ in range(a.permutations)])
        entry["reversal_score"]["residual_permutation_p"] = float((1 + (np.abs(null) >= abs(entry["reversal_score"]["spearman_residual"])).sum()) / (1 + a.permutations))
        entry["reversal_score"]["residual_null_sd"] = float(null.std())
        res[metric] = entry

    pairs.drop(columns="pair").to_csv(out / "pair_scores.tsv", sep="\t", index=False)
    (out / "summary.json").write_text(json.dumps(res, indent=2, default=float))
    print(f"pairs: {res['n_pairs']} ({res['n_drugs']} drugs, {res['n_combination_rows']} measured rows)")
    for metric in SYNERGY:
        e = res[metric]
        print(f"{metric:15s} reversal: rho_raw {e['reversal_score']['spearman_raw']:+.3f} "
              f"rho_resid {e['reversal_score']['spearman_residual']:+.3f} "
              f"(p={e['reversal_score']['residual_permutation_p']:.4f}) "
              f"AUC {e['reversal_score']['auc_synergy_vs_antagonism']:.3f} "
              f"[n+={e['reversal_score']['n_synergy']} n-={e['reversal_score']['n_antagonism']}]  "
              f"|cos| baseline: rho_resid {e['abs_cos_signatures']['spearman_residual']:+.3f} "
              f"AUC {e['abs_cos_signatures']['auc_synergy_vs_antagonism']:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
