"""LINCS Phase I consensus signatures via the CLUE API, for the drugs Phase II lacks.

Why this exists: Phase II (GSE70138) carries none of erastin, the platinums or most
vinca alkaloids, so this project built GEO consensus signatures instead -- and that
failed. No drug reached the 0.30 cross-series agreement the validation set for
"medium", and vincristine's two series *anti-correlate* (docs/drug_consensus.md).
Phase I has far more cell lines per drug (erastin 51, idarubicin 51, pemetrexed 50),
but its Level 5 matrix is 19.9 GB gzipped and cannot be partially fetched, because a
.gctx is HDF5 and needs random access while a gzip stream cannot be seeked.

The CLUE API serves the same signatures a few kB at a time. What it does NOT serve is
the gene-level z-vector: `/api/sigs` returns each signature's top 100 up and top 100
down genes (`up100_bing` / `dn100_bing`, Entrez ids in BING space). Probed and
confirmed -- there is no dataspace or l1000 endpoint on this key.

So a signature here is a *set*, not a vector, and the consensus is built accordingly:

  1. Within a cell line, a gene scores (fraction of that line's signatures putting it
     up) - (fraction putting it down). This averages over dose and timepoint.
  2. Across cell lines, the per-line scores are averaged. A gene reaching +1 was in
     the top 100 up in every signature of every cell line; a gene near 0 was
     idiosyncratic to one context.

Step 2 is the whole point. A LINCS consensus is trustworthy in this project because it
medians over cell lines, and the same drug in two cell lines otherwise agrees at only
0.18 (docs/npi_drug_retrieval.md). The output records `n_cell_lines` so the
reliability band follows from the data rather than from the source's reputation.

**These signatures are sparse and are NOT on the same scale as the dense Phase II
vectors.** Cosines against them cannot be compared with cosines against the Phase II
panel -- exactly the cross-panel error corrected in docs/drug_consensus.md. Compare
within this panel only.

The API key is read from CLUE_API_KEY and is never written to disk.

    export CLUE_API_KEY=...
    python scripts/ingest_clue_phase1.py
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))
from npi_pharma.model import DRUG, Signature  # noqa: E402
from npi_pharma.store import save_signatures  # noqa: E402

API = "https://api.clue.io/api/sigs"
FIELDS = ["sig_id", "pert_iname", "cell_id", "pert_idose", "pert_itime",
          "up100_bing", "dn100_bing", "distil_cc_q75"]

# The drugs Phase II lacks that Phase I actually carries, from
# scripts/scope_lincs_phase1.py. carboplatin, oxaliplatin, bleomycin, carmustine,
# lomustine and triapine are absent from Phase I and cannot be fetched at all.
DRUGS = ["erastin", "idarubicin", "pemetrexed", "vinblastine", "vincristine",
         "vinorelbine", "topotecan", "dactinomycin", "methotrexate",
         "cyclophosphamide", "cisplatin", "fluorouracil", "melphalan", "azacitidine"]


def fetch_sigs(drug: str, key: str, page: int = 1000, retries: int = 4) -> list[dict]:
    """Every trt_cp signature for one compound, paged, with the fields above."""
    out, skip = [], 0
    while True:
        flt = json.dumps({"where": {"pert_iname": drug, "pert_type": "trt_cp"},
                          "fields": FIELDS, "limit": page, "skip": skip})
        url = f"{API}?{urllib.parse.urlencode({'filter': flt})}"
        req = urllib.request.Request(url, headers={"user_key": key,
                                                   "Accept": "application/json"})
        for attempt in range(retries):
            try:
                with urllib.request.urlopen(req, timeout=180) as h:
                    got = json.loads(h.read().decode())
                break
            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as e:
                if attempt == retries - 1:
                    raise
                time.sleep(2 ** attempt)
                print(f"    retry {attempt + 1} for {drug}: {e}", file=sys.stderr)
        out += got
        if len(got) < page:
            return out
        skip += page


def consensus_from_sets(sigs: list[dict], symbol: dict[str, str]
                        ) -> tuple[pd.Series, dict]:
    """Cross-cell-line consensus over top-100 up/down sets, and what went into it.

    Per cell line a gene scores (up fraction) - (down fraction) over that line's
    signatures, so dose and timepoint are averaged within the line first and one
    heavily profiled line cannot outvote the rest. The per-line scores are then
    averaged across lines.
    """
    by_cell: dict[str, list[dict]] = defaultdict(list)
    for s in sigs:
        if s.get("up100_bing") and s.get("dn100_bing"):
            by_cell[s["cell_id"]].append(s)
    if not by_cell:
        return pd.Series(dtype=float), {"n_cell_lines": 0, "n_signatures": 0}

    per_cell = []
    for cell, group in by_cell.items():
        tally: dict[str, float] = defaultdict(float)
        for s in group:
            for g in s["up100_bing"]:
                tally[str(g)] += 1.0
            for g in s["dn100_bing"]:
                tally[str(g)] -= 1.0
        per_cell.append(pd.Series({k: v / len(group) for k, v in tally.items()}))

    M = pd.concat(per_cell, axis=1).fillna(0.0)
    score = M.mean(axis=1)
    score.index = [symbol.get(i, i) for i in score.index]
    score = score.groupby(level=0).mean()
    score = score[score.index.notna() & (score.index != "")]

    # How reproducible is the drug's own direction across cell lines? This is the
    # analogue of cross-series agreement for the GEO consensuses, computed the same
    # way (median pairwise cosine), so the two panels' reliability is comparable.
    Z = M.apply(lambda c: (c - c.mean()) / c.std(ddof=0) if c.std(ddof=0) > 0 else c * 0)
    C = (Z.T @ Z) / len(Z)
    iu = np.triu_indices_from(C, 1)
    pairs = C.to_numpy()[iu]
    ok = pairs[np.isfinite(pairs)]
    return score, {
        "n_cell_lines": int(M.shape[1]),
        "n_signatures": int(sum(len(g) for g in by_cell.values())),
        "cell_lines": sorted(by_cell),
        "median_cross_cell_line_cosine": float(np.median(ok)) if len(ok) else None,
        "n_genes": int(len(score)),
    }


def split_half_over_cells(sigs: list[dict], symbol: dict[str, str],
                          n_splits: int = 20, seed: int = 0) -> dict:
    """Does the consensus reproduce when built from two disjoint halves of the lines?

    This is the reliability that matters for a consensus, and it is not the pairwise
    agreement between cell lines: averaging over n units shrinks noise as 1/sqrt(n),
    so a consensus over 51 weakly-agreeing lines can still be highly reproducible.
    """
    cells = sorted({s["cell_id"] for s in sigs if s.get("up100_bing")})
    if len(cells) < 4:
        return {"n_units": len(cells), "split_half": None,
                "note": "needs >= 4 cell lines to split into two halves of >= 2"}
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_splits):
        perm = list(cells)
        rng.shuffle(perm)
        h = len(perm) // 2
        left, _ = consensus_from_sets([s for s in sigs if s["cell_id"] in set(perm[:h])],
                                      symbol)
        right, _ = consensus_from_sets(
            [s for s in sigs if s["cell_id"] in set(perm[h:2 * h])], symbol)
        g = left.index.intersection(right.index)
        if len(g) < 200 or left[g].std(ddof=0) == 0 or right[g].std(ddof=0) == 0:
            continue
        x = (left[g] - left[g].mean()) / left[g].std(ddof=0)
        y = (right[g] - right[g].mean()) / right[g].std(ddof=0)
        vals.append(float((x * y).mean()))
    if not vals:
        return {"n_units": len(cells), "split_half": None, "note": "no usable split"}
    return {"n_units": len(cells), "split_half": float(np.mean(vals)),
            "split_half_sd": float(np.std(vals)), "n_splits_used": len(vals)}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--drugs", nargs="*", default=DRUGS)
    ap.add_argument("--gene-info",
                    default="data/raw/GSE92742/GSE92742_Broad_LINCS_gene_info.txt.gz")
    ap.add_argument("--min-cell-lines", type=int, default=2,
                    help="a single cell line is not a consensus")
    ap.add_argument("--out", default="data/processed/signatures/lincs_phase1_clue.parquet")
    ap.add_argument("--report", default="out/lincs_phase1_clue")
    ap.add_argument("--split-halves", type=int, default=20,
                    help="random cell-line splits per drug for the reproducibility test; "
                         "0 skips it")
    a = ap.parse_args(argv)

    key = os.environ.get("CLUE_API_KEY", "").strip()
    if not key:
        print("CLUE_API_KEY is not set. Export it; it is never read from a file here.",
              file=sys.stderr)
        return 2
    rep = Path(a.report)
    rep.mkdir(parents=True, exist_ok=True)

    gi = pd.read_csv(a.gene_info, sep="\t", dtype=str)
    symbol = dict(zip(gi["pr_gene_id"], gi["pr_gene_symbol"]))

    sigs_out, log, splithalf = [], [], {}
    for drug in a.drugs:
        try:
            raw = fetch_sigs(drug, key)
        except Exception as e:
            log.append({"drug": drug, "outcome": f"error: {type(e).__name__}: {e}"})
            print(f"  {drug}: {type(e).__name__}", file=sys.stderr)
            continue
        score, qc = consensus_from_sets(raw, symbol)
        if qc["n_cell_lines"] < a.min_cell_lines:
            log.append({"drug": drug, "outcome": f"only {qc['n_cell_lines']} cell line(s)",
                        **{k: v for k, v in qc.items() if k != "cell_lines"}})
            print(f"  {drug}: only {qc['n_cell_lines']} cell line(s), skipped",
                  file=sys.stderr)
            continue
        sigs_out.append(Signature(
            sig_id=drug, kind=DRUG, genes=list(score.index),
            z=score.to_numpy(dtype=float),
            provenance="lincs_phase1_clue_api_top100_sets", modality="compound",
            species="human", sample_size=qc["n_cell_lines"],
            contrast="cross-cell-line consensus of top-100 up/down sets "
                     "(score = up fraction - down fraction, averaged over cell lines)",
            quality_flag="ok" if qc["n_cell_lines"] >= 7 else "few_cell_lines",
            source_accessions=["GSE92742"],
            meta={"drug_id": drug, "source": "CLUE API /api/sigs",
                  "sparse_not_comparable_to_dense_panels": True,
                  **{k: v for k, v in qc.items() if k != "cell_lines"},
                  "cell_lines": qc["cell_lines"]}))
        # Split-half reproducibility of the consensus itself, which is the statistic
        # that matters -- pairwise cross-cell-line agreement measures whether two lines
        # agree, not whether their average reproduces, and understates it badly (see
        # scripts/validate_consensus_splithalf.py). Computed here because it needs the
        # per-cell-line units, which are not kept after the consensus is taken.
        if a.split_halves:
            splithalf[drug] = split_half_over_cells(raw, symbol, a.split_halves)
        log.append({"drug": drug, "outcome": "built",
                    **{k: v for k, v in qc.items() if k != "cell_lines"},
                    "split_half": (splithalf.get(drug) or {}).get("split_half")})
        sh = (splithalf.get(drug) or {}).get("split_half")
        print(f"  {drug}: {qc['n_cell_lines']} cell lines, {qc['n_signatures']} signatures, "
              f"{qc['n_genes']} genes, pairwise {qc['median_cross_cell_line_cosine']:.3f}, "
              f"split-half {sh:+.3f}" if sh is not None else
              f"  {drug}: {qc['n_cell_lines']} cell lines, split-half n/a", file=sys.stderr)

    t = pd.DataFrame(log)
    t.to_csv(rep / "log.tsv", sep="\t", index=False)
    if splithalf:
        (rep / "splithalf.json").write_text(json.dumps(splithalf, indent=2, default=float))
    if sigs_out:
        save_signatures(sigs_out, a.out)
    (rep / "summary.json").write_text(json.dumps({
        "n_signatures": len(sigs_out),
        "scale_warning":
            "These are sparse set-derived signatures. Cosines against them are NOT on "
            "the same scale as cosines against the dense LINCS Phase II panel; compare "
            "within this panel only.",
    }, indent=2))
    print(f"\n{len(sigs_out)} Phase I consensus signatures -> {a.out}")
    built = t[t["outcome"] == "built"] if len(t) else t
    if len(built):
        print(built[["drug", "n_cell_lines", "n_signatures", "n_genes",
                     "median_cross_cell_line_cosine"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
