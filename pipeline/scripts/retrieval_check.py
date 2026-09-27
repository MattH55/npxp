"""Regression check: does the LCD signature retrieve the paper's L1000CDS2 mimics?

J Transl Med 2025 (doi:10.1186/s12967-025-07424-z) queried L1000CDS2 with the
GSE95640 LCD DEGs and reported the top 50 mimic signatures
(configs/retrieval/GSE95640_L1000CDS2_top50.tsv). This script scores every
trt_cp signature in a LINCS Level 5 GCTX by cosine similarity to our LCD
signature (cos > 0 = mimics LCD), then aggregates per compound (BRD id):

  median   median cosine over the compound's signatures (like ingest-lincs consensus)
  best     max cosine over its signatures (closer to L1000CDS2, which ranks
           single cell-line signatures)

It then asks where the paper's hits that exist in this release rank among all
compounds (AUC; permutation p against random compound sets of the same size).
The paper used Phase I (GSE92742) signatures and its own DEG list, so agreement
is expected to be partial.

    python scripts/retrieval_check.py --npis data/processed/signatures/npis.parquet \
        --lincs-dir data/raw/GSE70138 --out out/retrieval
"""

from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from npi_pharma.genes import normalize_symbols
from npi_pharma.ingest.lincs import GCTX_COL_IDS, GCTX_MATRIX, GCTX_ROW_IDS, _decode, load_gene_map
from npi_pharma.store import load_signatures


def _one(pattern: str) -> str:
    hits = sorted(glob.glob(pattern))
    if not hits:
        raise SystemExit(f"no file matches {pattern}")
    return hits[-1]


def auc_of(scores: pd.Series, positives: set[str]) -> float:
    r = scores.rank()
    pos = r[r.index.isin(positives)]
    n1, n0 = len(pos), len(r) - len(pos)
    return float((pos.sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--npis", default="data/processed/signatures/npis.parquet")
    ap.add_argument("--npi", default="LCD_adipose_GSE95640")
    ap.add_argument("--lincs-dir", default="data/raw/GSE70138")
    ap.add_argument("--hits", default="configs/retrieval/GSE95640_L1000CDS2_top50.tsv")
    ap.add_argument("--gene-space", choices=["landmark", "bing", "all"], default="bing")
    ap.add_argument("--block", type=int, default=4000)
    ap.add_argument("--permutations", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="out/retrieval")
    a = ap.parse_args(argv)

    npi = next(s for s in load_signatures(a.npis) if s.sig_id == a.npi).as_series()
    d = Path(a.lincs_dir)
    gctx = _one(str(d / "*Level5_COMPZ*.gctx"))
    si = pd.read_csv(_one(str(d / "*sig_info*.txt")), sep="\t", dtype=str)
    gmap = load_gene_map(_one(str(d / "*gene_info*.txt")), a.gene_space)

    with h5py.File(gctx, "r") as f:
        rows = _decode(f[GCTX_ROW_IDS][:])
        cols = _decode(f[GCTX_COL_IDS][:])
        sym = normalize_symbols([gmap.get(r, "") for r in rows])
        keep = [i for i, s in enumerate(sym) if s and s in npi.index]
        seen, uniq = set(), []
        for i in keep:  # one row per symbol
            if sym[i] not in seen:
                seen.add(sym[i])
                uniq.append(i)
        v = npi.loc[[sym[i] for i in uniq]].to_numpy(float)
        v = v / np.linalg.norm(v)
        si = si[si["pert_type"] == "trt_cp"].set_index("sig_id")
        pos = {c: i for i, c in enumerate(cols)}
        idx = sorted(pos[s] for s in si.index if s in pos)
        cos = np.empty(len(idx))
        dset = f[GCTX_MATRIX]
        for start in range(0, len(idx), a.block):
            chunk = idx[start:start + a.block]
            m = dset[chunk[0]:chunk[-1] + 1, :][np.array(chunk) - chunk[0]][:, uniq].astype(float)
            cos[start:start + len(chunk)] = (m @ v) / np.maximum(np.linalg.norm(m, axis=1), 1e-12)
            print(f"  {start + len(chunk)}/{len(idx)} signatures", file=sys.stderr)
    per_sig = pd.DataFrame({"cos": cos}, index=[cols[i] for i in idx]).join(si[["pert_id", "pert_iname", "cell_id"]])
    per_sig["brd"] = per_sig["pert_id"].str[:13]

    comp = per_sig.groupby("brd").agg(name=("pert_iname", "first"), n_sigs=("cos", "size"),
                                      median=("cos", "median"), best=("cos", "max"))
    hits = pd.read_csv(a.hits, sep="\t", comment="#", dtype=str)
    paper = set(hits["pert_id"]) & set(comp.index)
    rng = np.random.default_rng(a.seed)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    comp.sort_values("best", ascending=False).to_csv(out / "compound_mimicry.tsv", sep="\t")
    print(f"LCD vs {len(comp)} {Path(gctx).name.split('_')[0]} compounds on {len(uniq)} {a.gene_space} genes; "
          f"paper hits present: {len(paper)}/{hits['pert_id'].nunique()}")
    for agg in ("median", "best"):
        s = comp[agg]
        auc = auc_of(s, paper)
        null = np.array([auc_of(s, set(rng.choice(comp.index, len(paper), replace=False)))
                         for _ in range(a.permutations)])
        p = (1 + (null >= auc).sum()) / (1 + a.permutations)
        pct = s.rank(pct=True)
        print(f"  {agg:6s} AUC {auc:.3f} (permutation p = {p:.4f}); "
              f"paper hits in our top 10%: {(pct[list(paper)] >= 0.9).sum()}/{len(paper)}")
    tbl = comp.loc[sorted(paper)].assign(best_pct=comp["best"].rank(pct=True), median_pct=comp["median"].rank(pct=True))
    print(tbl.sort_values("best_pct", ascending=False).round(3).to_string())
    print("our top 15 by best-signature mimicry:")
    print(comp.sort_values("best", ascending=False).head(15).round(3).to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
