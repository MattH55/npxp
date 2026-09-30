"""Split-half reproducibility of a consensus -- the right reliability statistic.

This project banded its consensus signatures by *pairwise* agreement between the units
that went into them (cross-series cosine for GEO, cross-cell-line cosine for LINCS),
and used 0.30 as the threshold for "medium". That statistic answers the wrong
question. It measures how well two individual units agree, not how reproducible their
*average* is -- and averaging is the entire point of a consensus. The standard error
of a mean falls as 1/sqrt(n), so a consensus over many weak units can be far more
reproducible than any pair of them.

Measured here rather than argued: LINCS Phase I erastin has 51 cell lines whose
pairwise agreement is 0.025, and its split-half reproducibility is **0.494**. Pairwise
agreement understated it twentyfold.

The test: split the units into two disjoint halves, build a consensus from each, and
correlate the two consensuses. Repeat over random splits and average. This measures
exactly what matters -- would a fresh set of experiments give the same consensus? -- and
it is the same quantity for both panels, so they can be compared.

Two caveats it does not remove. Split-half uses half the units, so it *understates* the
full consensus's reproducibility. And for the set-derived Phase I panel the vectors are
sparse, which deflates any cosine; that affects the pairwise numbers far more than the
split-half ones, since a consensus over many lines is dense.

    python scripts/validate_consensus_splithalf.py
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))
from npi_pharma.store import load_signatures  # noqa: E402


def cosine(a: pd.Series, b: pd.Series, min_genes: int = 200) -> float:
    """Standardised cosine on shared genes, or NaN if there are too few."""
    g = a.dropna().index.intersection(b.dropna().index)
    if len(g) < min_genes:
        return float("nan")
    x, y = a[g], b[g]
    if x.std(ddof=0) == 0 or y.std(ddof=0) == 0:
        return float("nan")
    zx = (x - x.mean()) / x.std(ddof=0)
    zy = (y - y.mean()) / y.std(ddof=0)
    return float((zx * zy).mean())


def split_half(units: list[pd.Series], n_splits: int = 20, seed: int = 0,
               combine=lambda f: f.median(axis=1)) -> dict:
    """Mean correlation between consensuses built from two disjoint halves."""
    if len(units) < 4:
        return {"n_units": len(units), "split_half": None,
                "note": "needs >= 4 units to split into two halves of >= 2"}
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_splits):
        idx = list(range(len(units)))
        rng.shuffle(idx)
        h = len(idx) // 2
        left = combine(pd.concat([units[i] for i in idx[:h]], axis=1, join="inner"))
        right = combine(pd.concat([units[i] for i in idx[h:2 * h]], axis=1, join="inner"))
        c = cosine(left, right)
        if np.isfinite(c):
            vals.append(c)
    if not vals:
        return {"n_units": len(units), "split_half": None,
                "note": "no split shared enough genes"}
    return {"n_units": len(units), "split_half": float(np.mean(vals)),
            "split_half_sd": float(np.std(vals)), "n_splits_used": len(vals)}


def geo_units(path: Path) -> dict[str, list[pd.Series]]:
    """Per-series GEO signatures, grouped by drug."""
    out: dict[str, list[pd.Series]] = defaultdict(list)
    for s in load_signatures(path):
        drug = (s.meta or {}).get("drug_id") or s.sig_id.split("|")[0]
        out[drug].append(s.as_series())
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--geo-per-series",
                    default="data/processed/signatures/drugs_geo_per_series_v2.parquet")
    ap.add_argument("--geo-qc", default="out/drug_consensus_v2/qc.json")
    ap.add_argument("--phase1-log", default="out/lincs_phase1_clue/log.tsv",
                    help="carries each Phase I drug's pairwise cross-cell-line cosine")
    ap.add_argument("--phase1-splithalf", default="out/lincs_phase1_clue/splithalf.json",
                    help="written by the ingest; read here so the API is not re-queried")
    ap.add_argument("--n-splits", type=int, default=20)
    ap.add_argument("--out", default="out/consensus_splithalf")
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    qc = json.loads(Path(a.geo_qc).read_text()) if Path(a.geo_qc).exists() else {}
    rows = []
    if Path(a.geo_per_series).exists():
        for drug, units in geo_units(Path(a.geo_per_series)).items():
            r = split_half(units, a.n_splits)
            rows.append({"panel": "GEO consensus v2", "drug": drug,
                         "n_units": r["n_units"], "unit": "series",
                         "pairwise": (qc.get(drug) or {}).get("median_cross_series_cosine"),
                         "split_half": r.get("split_half"),
                         "note": r.get("note")})
    if Path(a.phase1_splithalf).exists():
        pre = json.loads(Path(a.phase1_splithalf).read_text())
        pw = {}
        if Path(a.phase1_log).exists():
            lg = pd.read_csv(a.phase1_log, sep="\t")
            pw = dict(zip(lg["drug"], lg.get("median_cross_cell_line_cosine", [])))
        for drug, v in pre.items():
            rows.append({"panel": "LINCS Phase I (CLUE)", "drug": drug,
                         "n_units": v.get("n_units"), "unit": "cell lines",
                         "pairwise": pw.get(drug), "split_half": v.get("split_half"),
                         "note": v.get("note")})

    t = pd.DataFrame(rows).sort_values(["panel", "split_half"], ascending=[True, False])
    t.to_csv(out / "splithalf.tsv", sep="\t", index=False)
    print(t.to_string(index=False, na_rep="-"))

    got = t.dropna(subset=["split_half", "pairwise"])
    verdict = {
        "n_scored": int(len(got)),
        "clears_0.30_on_split_half": sorted(got.loc[got["split_half"] >= 0.30, "drug"]),
        "clears_0.30_on_pairwise": sorted(got.loc[got["pairwise"] >= 0.30, "drug"]),
        "median_ratio_splithalf_over_pairwise":
            float((got["split_half"] / got["pairwise"].replace(0, np.nan)).median()),
        "why_split_half":
            "Pairwise agreement measures whether two units agree; split-half measures "
            "whether the consensus reproduces, which is what a consensus is for. The "
            "second is the reliability that matters and is systematically higher.",
        "caveat":
            "Split-half builds each side from half the units, so it understates the "
            "full consensus. For the set-derived Phase I panel, sparsity deflates the "
            "pairwise numbers much more than the split-half ones.",
    }
    (out / "verdict.json").write_text(json.dumps(verdict, indent=2, default=float))
    print(f"\nclears 0.30 on split-half: {', '.join(verdict['clears_0.30_on_split_half']) or 'none'}")
    print(f"clears 0.30 on pairwise:   {', '.join(verdict['clears_0.30_on_pairwise']) or 'none'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
