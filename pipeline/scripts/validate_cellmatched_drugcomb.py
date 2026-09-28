"""Strongest available test: signatures measured in the SAME cell line as the synergy.

The pair-level and cell-wise tests both use LINCS consensus signatures pooled
across cell lines, which leaves one loophole: maybe signature composition works,
but only when the signature comes from the same cell line as the measured
combination. Nine cell lines are in both DrugComb and LINCS Phase II, giving
~3.4k (pair, cell line) observations where both drugs were profiled in that
exact line.

For each such observation the npi_pharma scores are computed with cell-matched
signatures, and tested against measured synergy after removing pair and
cell-line main effects.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from npi_pharma.cancer.programs import center_within, load_expression
from npi_pharma.genes import normalize_symbols
from npi_pharma.ingest.lincs import load_gene_map, read_gctx_columns

SYNERGY = ["synergy_bliss", "synergy_loewe", "synergy_zip", "synergy_hsa"]
NORM = lambda s: re.sub(r"[^a-z0-9]", "", str(s).lower())


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--drugcomb", required=True)
    ap.add_argument("--lincs-dir", default="data/raw/GSE70138")
    ap.add_argument("--expression", default="data/raw/depmap/OmicsExpressionProteinCodingGenesTPMLogp1.csv")
    ap.add_argument("--models", default="data/raw/depmap/Model.csv")
    ap.add_argument("--focus-top-k", type=int, default=500)
    ap.add_argument("--permutations", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="out/drugcomb_validation")
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    ld = Path(a.lincs_dir)

    d = pd.read_csv(a.drugcomb, usecols=["drug_row", "drug_col", "cell_line_name", *SYNERGY], low_memory=False)
    d = d[d["drug_col"].notna() & (d["drug_col"] != "NULL")]
    for c in SYNERGY:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.dropna(subset=["synergy_bliss"])

    si = pd.read_csv(next(ld.glob("*sig_info*.txt")), sep="\t",
                     usecols=["sig_id", "pert_iname", "pert_type", "cell_id"])
    si = si[si["pert_type"] == "trt_cp"]
    names = {NORM(n): n for n in si["pert_iname"].unique()}
    cells = {NORM(c): c for c in si["cell_id"].unique()}
    d["A"] = d["drug_row"].map(lambda x: names.get(NORM(x)))
    d["B"] = d["drug_col"].map(lambda x: names.get(NORM(x)))
    d["C"] = d["cell_line_name"].map(lambda x: cells.get(NORM(x)))
    d = d.dropna(subset=["A", "B", "C"])
    d = d[d["A"] != d["B"]]
    have = set(zip(si["pert_iname"], si["cell_id"]))
    d = d[[(r.A, r.C) in have and (r.B, r.C) in have for r in d.itertuples()]]
    d["pair"] = [" + ".join(sorted(p)) for p in zip(d["A"], d["B"])]
    d = d.groupby(["pair", "A", "B", "C"], as_index=False)[SYNERGY].median()
    need = set(zip(d["A"], d["C"])) | set(zip(d["B"], d["C"]))
    print(f"{len(d)} (pair, cell line) observations; {d['pair'].nunique()} pairs; "
          f"{d['C'].nunique()} cell lines; {len(need)} cell-matched signatures needed", file=sys.stderr)

    sel = si[[(r.pert_iname, r.cell_id) in need for r in si.itertuples()]]
    gmap = load_gene_map(next(ld.glob("*gene_info*.txt")), "landmark")
    block = read_gctx_columns(next(ld.glob("*Level5_COMPZ*.gctx")), sel["sig_id"].tolist())
    sym = normalize_symbols([gmap.get(i, "") for i in block.index])
    block = block[[s != "" for s in sym]]
    block.index = [s for s in sym if s != ""]
    key = sel.set_index("sig_id")
    cols = pd.MultiIndex.from_arrays([key.loc[block.columns, "pert_iname"], key.loc[block.columns, "cell_id"]])
    block.columns = cols
    sigs = block.T.groupby(level=[0, 1]).median().T            # genes x (drug, cell)
    sigs = sigs.groupby(level=0).median()   # one row per symbol
    print(f"{sigs.shape[0]} landmark genes x {sigs.shape[1]} cell-matched signatures", file=sys.stderr)

    models = pd.read_csv(a.models)
    lineage = models.set_index("ModelID")["OncotreeLineage"]
    mm = {NORM(r.StrippedCellLineName): r.ModelID for r in models.itertuples()
          if isinstance(r.StrippedCellLineName, str)}
    expr = load_expression(a.expression, list(sigs.index))
    genes = sigs.index.intersection(expr.columns)
    sigs = sigs.loc[genes]
    xc = center_within(expr[genes], lineage)
    P = ((xc - xc.mean()) / xc.std(ddof=0)).T
    d["ModelID"] = d["C"].map(lambda c: mm.get(NORM(c)))
    d = d[d["ModelID"].isin(P.columns)]
    print(f"{len(d)} observations with DepMap expression for the cell line", file=sys.stderr)

    unit = sigs / np.sqrt((sigs ** 2).sum())
    k = min(a.focus_top_k, len(genes))
    rows = []
    for r in d.itertuples():
        p = P[r.ModelID].to_numpy()
        foc = np.argsort(-np.abs(p))[:k]
        pf = p[foc]
        sa, sb = unit[(r.A, r.C)].to_numpy(), unit[(r.B, r.C)].to_numpy()
        rv = lambda x: float(-(pf @ x[foc]) / (np.linalg.norm(pf) * np.linalg.norm(x[foc])))
        r_a, r_b, r_c = rv(sa), rv(sb), rv(sa + sb)
        best = max(r_a, r_b)
        rows.append({"pair": r.pair, "cell": r.C, "mono_reversal": (r_a + r_b) / 2,
                     "complementarity_gain": r_c - best,
                     "complementarity": (r_c - best) / max(1 - best, 1e-6),
                     "reverse_combo": r_c, **{c: getattr(r, c) for c in SYNERGY}})
    t = pd.DataFrame(rows)
    t.to_csv(out / "cellmatched_scores.tsv", sep="\t", index=False)

    feats = ["complementarity", "complementarity_gain", "mono_reversal", "reverse_combo"]
    res = {"n_observations": len(t), "n_pairs": int(t["pair"].nunique()), "n_cell_lines": int(t["cell"].nunique()),
           "n_genes": int(len(genes)), "focus_top_k": a.focus_top_k,
           "cell_lines": t["cell"].value_counts().to_dict()}
    rng = np.random.default_rng(a.seed)
    for metric in SYNERGY:
        y = t[metric]
        resid = (y - t.groupby("pair")[metric].transform("mean")
                 - t.groupby("cell")[metric].transform("mean") + y.mean()).to_numpy()
        e = {f: {"spearman_raw": float(t[f].corr(y, method="spearman")),
                 "spearman_residual": float(t[f].corr(pd.Series(resid), method="spearman"))} for f in feats}
        obs = e["complementarity"]["spearman_residual"]
        null = np.array([float(pd.Series(rng.permutation(t["complementarity"].to_numpy())).corr(pd.Series(resid), method="spearman"))
                         for _ in range(a.permutations)])
        e["complementarity"]["residual_permutation_p"] = float((1 + (np.abs(null) >= abs(obs)).sum()) / (1 + a.permutations))
        e["complementarity"]["residual_null_sd"] = float(null.std())
        res[metric] = e
    (out / "cellmatched_summary.json").write_text(json.dumps(res, indent=2, default=float))
    print(f"observations: {res['n_observations']}  pairs: {res['n_pairs']}  cell lines: {res['n_cell_lines']}")
    for metric in SYNERGY:
        e = res[metric]
        print(f"{metric:15s} " + "  ".join(
            f"{f}: raw {e[f]['spearman_raw']:+.3f} resid {e[f]['spearman_residual']:+.3f}" for f in feats)
            + f"  (complementarity p={e['complementarity']['residual_permutation_p']:.4f}, "
              f"null sd {e['complementarity']['residual_null_sd']:.3f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
