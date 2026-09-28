"""Find GEO expression datasets where a drug was applied to human cells.

LINCS L1000 is not enough on its own: Phase II (GSE70138) contains only 5 of this
project's 13 curated drugs, and none of the platinums (cisplatin, carboplatin,
oxaliplatin), 5-fluorouracil, erastin, triapine, lomustine or cyclophosphamide.
Drug signatures for those have to come from ordinary GEO series.

This searches GEO DataSets for human expression series naming a drug, scores each
by how likely it is to be a usable treated-vs-control cell-line experiment, and
writes a ranked candidate table. It does not build signatures: whether a series
really has a drug arm and a matched control has to be read off its sample
annotations (`npi-pharma inspect-geo`), the same way the NPI catalog is curated.

    python scripts/find_drug_datasets.py --drugs cisplatin oxaliplatin erastin
    python scripts/find_drug_datasets.py --missing-from-lincs --max-per-drug 120
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

CURATED_DRUGS = ["cisplatin", "carboplatin", "oxaliplatin", "5-fluorouracil", "mitomycin C",
                 "doxorubicin", "paclitaxel", "temozolomide", "metformin", "erastin", "triapine",
                 "lomustine", "cyclophosphamide"]
# Present in LINCS Phase II (GSE70138), so lower priority for GEO curation.
IN_LINCS_PHASE2 = {"mitomycin C", "doxorubicin", "paclitaxel", "temozolomide", "metformin"}

CURATED_CELL_LINES = [
    "HeLa", "CaSki", "MCF-7", "MCF7", "T47D", "T-47D", "MDA-MB-231", "MDA-MB-468", "A549",
    "U-87", "U87", "T98G", "DU145", "DU-145", "U937", "U-937", "RKO", "HCT116", "HCT-116",
    "LoVo", "Hs578T", "HT29", "HT-29", "SW480", "SW620", "PC3", "PC-3", "A375", "HEPG2", "HepG2",
]

CELL_CONTEXT = re.compile(r"cell line|cell lines|in vitro|cultured|cells were treated", re.I)
TREATED = re.compile(r"treat|expos|incubat|versus control|vs\.? control|untreated|dose|concentration|"
                     r"\buM\b|µM|μM|\bnM\b|mg/ml|resistan", re.I)
PATIENT = re.compile(r"patient|biopsy|tumor tissue|tumour tissue|xenograft|blood sample|serum|"
                     r"clinical trial|cohort|responder", re.I)
RESISTANT_ONLY = re.compile(r"resistant (cell )?lines? (were|compared)|acquired resistance", re.I)


def get(url: str, timeout: int = 90, tries: int = 5) -> bytes:
    """GET with backoff. NCBI allows ~3 requests/second unauthenticated and answers
    HTTP 429 above that, so a burst of esummary calls must back off rather than fail."""
    delay = 1.0
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code not in (429, 500, 502, 503) or attempt == tries - 1:
                raise
        except (urllib.error.URLError, TimeoutError):
            if attempt == tries - 1:
                raise
        time.sleep(delay)
        delay *= 2
    raise RuntimeError("unreachable")


def esearch(term: str, retmax: int) -> list[str]:
    q = urllib.parse.urlencode({"db": "gds", "term": term, "retmax": retmax, "retmode": "json"})
    return json.loads(get(f"{EUTILS}/esearch.fcgi?{q}")).get("esearchresult", {}).get("idlist", [])


def esummary(uids: list[str]) -> list[dict]:
    q = urllib.parse.urlencode({"db": "gds", "id": ",".join(uids), "retmode": "json"})
    d = json.loads(get(f"{EUTILS}/esummary.fcgi?{q}")).get("result", {})
    return [d[u] for u in d.get("uids", [])]


PROCESSED = re.compile(r"\.(txt|tsv|csv|xlsx?)(\.gz)?$|fpkm|tpm|count|rma|normali[sz]ed|matrix", re.I)


def has_processed_supplement(gse: str, timeout: int = 30) -> str | None:
    """Names of a series' supplementary files that look like a processed expression
    table, or None. A series whose only supplement is *_RAW.tar needs raw array or
    fastq processing, which is out of scope here."""
    stem = re.match(r"([A-Z]+)(\d+)$", gse)
    if not stem:
        return None
    pre, dig = stem.groups()
    bucket = f"{pre}{dig[:-3]}nnn" if len(dig) > 3 else f"{pre}nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{bucket}/{gse}/suppl/"
    try:
        html = get(url, timeout).decode("utf-8", "replace")
    except Exception:
        return None
    names = {h for h in re.findall(r'href="([^"?/][^"]*)"', html) if not h.startswith("http")}
    good = sorted(n for n in names if PROCESSED.search(n) and not n.endswith("_RAW.tar")
                  and "filelist" not in n.lower())
    return ",".join(good[:3]) or None


def score(rec: dict, drug: str) -> tuple[int, list[str]]:
    """Heuristic priority. Higher is a better candidate for a clean drug signature."""
    text = f"{rec.get('title', '')} {rec.get('summary', '')}"
    n = int(rec.get("n_samples") or 0)
    s, why = 0, []
    if CELL_CONTEXT.search(text):
        s += 3
        why.append("cell-line context")
    if TREATED.search(text):
        s += 2
        why.append("treatment wording")
    if PATIENT.search(text):
        s -= 3
        why.append("patient/xenograft wording")
    if RESISTANT_ONLY.search(text):
        s -= 1
        why.append("resistance-line comparison")
    hits = sorted({c for c in CURATED_CELL_LINES if re.search(rf"\b{re.escape(c)}\b", text)})
    if hits:
        s += 2
        why.append("curated line: " + ",".join(hits[:4]))
    if 4 <= n <= 40:
        s += 2
        why.append(f"n={n}")
    elif 41 <= n <= 120:
        s += 1
        why.append(f"n={n}")
    elif n > 300:
        s -= 2
        why.append(f"large n={n}")
    if re.search(rf"\b{re.escape(drug)}\b", rec.get("title", ""), re.I):
        s += 2
        why.append("drug in title")
    if "high throughput sequencing" in (rec.get("gdstype") or "").lower():
        s += 1
        why.append("RNA-seq")
    return s, why


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--drugs", nargs="*", help="drug names to search (default: all curated)")
    ap.add_argument("--missing-from-lincs", action="store_true",
                    help="only the curated drugs absent from LINCS Phase II")
    ap.add_argument("--max-per-drug", type=int, default=80)
    ap.add_argument("--batch", type=int, default=100)
    ap.add_argument("--min-score", type=int, default=5)
    ap.add_argument("--check-supplement", action="store_true",
                   help="for the shortlist, check whether a processed expression table is published")
    ap.add_argument("--out", default="out/drug_dataset_candidates")
    a = ap.parse_args(argv)
    drugs = a.drugs or [d for d in CURATED_DRUGS if not a.missing_from_lincs or d not in IN_LINCS_PHASE2]

    rows = []
    for drug in drugs:
        term = (f'"{drug}" AND "Homo sapiens"[orgn] AND gse[etyp] AND '
                f'("expression profiling by array"[gdstyp] OR '
                f'"expression profiling by high throughput sequencing"[gdstyp])')
        uids = esearch(term, a.max_per_drug)
        print(f"{drug}: {len(uids)} series", file=sys.stderr)
        for i in range(0, len(uids), a.batch):
            for rec in esummary(uids[i:i + a.batch]):
                sc, why = score(rec, drug)
                rows.append({
                    "drug": drug, "accession": rec.get("accession"), "score": sc,
                    "n_samples": int(rec.get("n_samples") or 0),
                    "gpl": (rec.get("gpl") or "").split(";")[0],
                    "type": (rec.get("gdstype") or "")[:34],
                    "title": (rec.get("title") or "")[:130],
                    "reasons": "; ".join(why),
                    "summary": (rec.get("summary") or "")[:300],
                })
            time.sleep(0.5)
        time.sleep(0.5)

    t = pd.DataFrame(rows).drop_duplicates(subset=["drug", "accession"])
    t = t.sort_values(["drug", "score", "n_samples"], ascending=[True, False, True])
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t.to_csv(out / "all_candidates.tsv", sep="\t", index=False)
    top = t[t["score"] >= a.min_score].copy()
    if a.check_supplement:
        print(f"checking supplementary files for {len(top)} shortlisted series...", file=sys.stderr)
        supp = []
        for acc in top["accession"]:
            supp.append(has_processed_supplement(acc))
            time.sleep(0.4)
        top["processed_supplement"] = supp
        top = top.sort_values(["drug", "processed_supplement", "score"],
                              ascending=[True, False, False], na_position="last")
    top.to_csv(out / "shortlist.tsv", sep="\t", index=False)
    print(f"\n{len(t)} series scored; {len(top)} at score >= {a.min_score} -> {out}/shortlist.tsv")
    print(t.groupby("drug")["score"].agg(["count", "max"]).to_string())
    for drug, g in top.groupby("drug"):
        print(f"\n== {drug} ==")
        print(g.head(6)[["accession", "score", "n_samples", "gpl", "title"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
