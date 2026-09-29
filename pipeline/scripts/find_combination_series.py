"""Find GEO series with a factorial design, for MEASURED interaction signatures.

Everything else in this project infers interaction from single-agent signatures,
and that was measured not to work: signature composition carries no detectable
information about synergy across 50k drug-pair observations
(docs/validation_drugcomb.md). A series that profiles control, agent A, agent B
AND the A+B combination is different in kind — the interaction can be *computed*
rather than predicted:

    I(A, B) = (combo - control) - [(A - control) + (B - control)]

which is the transcriptional excess over additivity, the expression analogue of a
Bliss excess. Those series are the only route in this repository to interaction
data that is observed rather than assumed.

This scans GEO for such designs and reports the arms it finds. It builds nothing:
which value is the control, which are single agents and which is the combination
is read off GEO's sample annotations, and anything that does not resolve cleanly
is reported rather than guessed.

    python scripts/find_combination_series.py --max-candidates 80
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from build_drug_consensus import CONTROL, RESISTANCE_FIELD, RESISTANCE_VALUE, SKIP_FIELD  # noqa: E402
from find_drug_datasets import esearch, esummary  # noqa: E402

from npi_pharma.ingest.fetch import fetch_geo_series  # noqa: E402
from npi_pharma.ingest.geo import read_series_matrix  # noqa: E402

COMBO = re.compile(r"\+|\bplus\b|combination|combined|combi\b|combo|co-?treat|dual|"
                   r"\bwith\b|\band\b|/", re.I)

QUERIES = [
    '("drug combination" OR "combination treatment" OR "combination therapy") AND '
    '("cell line" OR "cell lines")',
    '(synergy OR synergistic OR synergism) AND ("cell line" OR "cell lines")',
    '(combination OR cotreatment OR "co-treatment") AND (chemotherapy OR cisplatin OR '
    'doxorubicin OR paclitaxel OR "5-fluorouracil" OR oxaliplatin)',
]
BASE = ('"Homo sapiens"[orgn] AND gse[etyp] AND ("expression profiling by array"[gdstyp] '
        'OR "expression profiling by high throughput sequencing"[gdstyp])')


def classify_arms(values: list[str]) -> dict[str, list[str]] | None:
    """Split a field's values into control / single agents / combination arms."""
    ctrl = [v for v in values if CONTROL.search(v) and not COMBO.search(v)]
    combo = [v for v in values if COMBO.search(v) and not CONTROL.search(v)]
    singles = [v for v in values if v not in ctrl and v not in combo]
    if len(ctrl) >= 1 and len(singles) >= 2 and len(combo) >= 1:
        return {"control": ctrl, "singles": singles, "combination": combo}
    return None


def detect_factorial(samples: pd.DataFrame, min_per_arm: int = 2) -> dict | None:
    """The best factorial field in a series' sample annotations, if any."""
    best = None
    for field in samples.columns:
        if field in SKIP_FIELD or RESISTANCE_FIELD.search(field):
            continue
        vals = samples[field].astype(str)
        uniq = list(vals.unique())
        if not (3 < len(uniq) <= 14) and len(uniq) < 4:
            continue
        if any(RESISTANCE_VALUE.search(v) for v in uniq):
            continue
        arms = classify_arms(uniq)
        if not arms:
            continue
        counts = vals.value_counts()
        if any(counts.get(v, 0) < min_per_arm
               for grp in arms.values() for v in grp):
            continue
        n = int(sum(counts.get(v, 0) for grp in arms.values() for v in grp))
        cand = {"field": field, **arms, "n_samples_used": n,
                "counts": {v: int(counts.get(v, 0)) for grp in arms.values() for v in grp}}
        if best is None or cand["n_samples_used"] > best["n_samples_used"]:
            best = cand
    return best


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--max-candidates", type=int, default=80, help="series per query to screen")
    ap.add_argument("--max-fetch", type=int, default=120, help="series matrices to actually fetch")
    ap.add_argument("--raw-dir", default="data/raw")
    ap.add_argument("--out", default="out/combination_series")
    a = ap.parse_args(argv)
    raw, out = Path(a.raw_dir), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    uids: list[str] = []
    for q in QUERIES:
        got = esearch(f"({q}) AND {BASE}", a.max_candidates)
        print(f"  {len(got)} series for: {q[:60]}...", file=sys.stderr)
        uids += got
        time.sleep(0.5)
    uids = list(dict.fromkeys(uids))
    recs = []
    for i in range(0, len(uids), 100):
        recs += esummary(uids[i:i + 100])
        time.sleep(0.5)
    # screen on the summary text before spending a fetch on each
    recs = [r for r in recs if COMBO.search(f"{r.get('title','')} {r.get('summary','')}")
            and 4 <= int(r.get("n_samples") or 0) <= 80]
    recs.sort(key=lambda r: int(r.get("n_samples") or 0))
    print(f"{len(uids)} unique series, {len(recs)} worth fetching", file=sys.stderr)

    rows = []
    for r in recs[: a.max_fetch]:
        acc = r.get("accession")
        try:
            fetch_geo_series(acc, raw, log=lambda *a_, **k: None)
            mats = sorted((raw / acc).glob("*series_matrix.txt.gz"))
            if not mats:
                continue
            _, samples = read_series_matrix(mats[0])
        except Exception as e:
            rows.append({"accession": acc, "outcome": f"error: {type(e).__name__}", "title": r.get("title", "")[:110]})
            continue
        hit = detect_factorial(samples)
        rows.append({
            "accession": acc, "n_samples": int(r.get("n_samples") or 0),
            "outcome": "factorial" if hit else "no factorial design detected",
            "field": hit["field"] if hit else None,
            "control": "; ".join(hit["control"])[:80] if hit else None,
            "singles": "; ".join(hit["singles"])[:120] if hit else None,
            "combination": "; ".join(hit["combination"])[:120] if hit else None,
            "counts": json.dumps(hit["counts"]) if hit else None,
            "title": r.get("title", "")[:110],
        })
        if hit:
            print(f"  {acc}: {len(hit['singles'])} singles + {len(hit['combination'])} combo "
                  f"on '{hit['field']}'", file=sys.stderr)
        time.sleep(0.2)

    t = pd.DataFrame(rows)
    t.to_csv(out / "candidates.tsv", sep="\t", index=False)
    hits = t[t["outcome"] == "factorial"]
    hits.to_csv(out / "factorial_series.tsv", sep="\t", index=False)
    print(f"\n{len(t)} series screened; {len(hits)} have a factorial design "
          f"-> {out}/factorial_series.tsv")
    if len(hits):
        with pd.option_context("display.width", 220, "display.max_colwidth", 46):
            print(hits[["accession", "n_samples", "field", "singles", "combination"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
