"""Search PMC open-access full text for NPI single-agent survival numbers.

Feeds `configs/npi_monotherapy.yaml`. The curation bar there is a published
surviving fraction or viability loss for the NPI ALONE in a named cell line, so
this script finds and prints *candidate sentences* with their citation; a human
(or Claude) then reads each one and decides. It never writes the YAML itself,
because deciding whether a sentence really states single-agent survival for a
named line is a judgement call, and a regex that guessed would put fabricated
numbers into the evidence base.

Uses NCBI E-utilities (esearch/efetch on the pmc database), which serve the
open-access subset without the reCAPTCHA gate on the PMC website.

    python scripts/find_npi_monotherapy.py --npi hyperthermia --max-papers 60
    python scripts/find_npi_monotherapy.py --npi glucose_restriction
    python scripts/find_npi_monotherapy.py --query "cystine deprivation surviving fraction" --cells MDAMB231
"""

from __future__ import annotations

import argparse
import html
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Cell lines this project curates, with the spellings papers actually use.
CELL_ALIASES = {
    "HELA_CERVIX": ["HeLa"],
    "CASKI_CERVIX": ["CaSki", "Ca Ski"],
    "MCF7_BREAST": ["MCF-7", "MCF7"],
    "T47D_BREAST": ["T47D", "T-47D"],
    "MDAMB231_BREAST": ["MDA-MB-231", "MDAMB231"],
    "MDAMB468_BREAST": ["MDA-MB-468", "MDAMB468"],
    "A549_LUNG": ["A549"],
    "U87MG_CENTRAL_NERVOUS_SYSTEM": ["U-87", "U87"],
    "T98G_CENTRAL_NERVOUS_SYSTEM": ["T98G", "T98-G"],
    "DU145_PROSTATE": ["DU145", "DU-145"],
    "U937_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE": ["U937", "U-937"],
    "RKO_LARGE_INTESTINE": ["RKO"],
    "HCT116_LARGE_INTESTINE": ["HCT116", "HCT-116"],
    "LOVO_LARGE_INTESTINE": ["LoVo", "LOVO"],
    "HS578T_BREAST": ["Hs578T", "Hs 578T"],
}

QUERIES = {
    "hyperthermia": [
        '(hyperthermia OR "heat shock" OR "thermal dose") AND ("surviving fraction" OR "clonogenic survival") '
        'AND ("cell line" OR "cell lines") AND open access[filter]',
        '(hyperthermia AND ("43 degrees" OR "43 C" OR "42 degrees")) AND (clonogenic OR "surviving fraction" '
        'OR viability) AND open access[filter]',
        '("thermal enhancement" OR thermosensitivity OR thermotolerance) AND ("surviving fraction" OR clonogenic) '
        'AND open access[filter]',
    ],
    "glucose_restriction": [
        '("glucose restriction" OR "low glucose" OR "glucose deprivation") AND (viability OR "surviving fraction" '
        'OR clonogenic) AND ("cancer cell" OR "cell line") AND open access[filter]',
    ],
    "amino_acid_restriction": [
        '("cystine deprivation" OR "cysteine deprivation" OR "methionine restriction" OR "methionine deprivation") '
        'AND (viability OR "cell death" OR clonogenic) AND ("cell line" OR "cancer cells") AND open access[filter]',
    ],
    "serum_starvation": [
        '("serum starvation" OR "serum-free" OR "serum deprivation") AND (viability OR clonogenic OR "surviving '
        'fraction") AND ("cancer cell line" OR "colon cancer") AND open access[filter]',
    ],
    "acidosis_hypoxia": [
        '(acidosis OR "low pH" OR hypoxia) AND (clonogenic OR "surviving fraction") AND (glioma OR glioblastoma) '
        'AND open access[filter]',
    ],
}

# A candidate sentence needs a number that could be a survival/viability value
NUMBER = re.compile(r"\d+(?:\.\d+)?\s?%|\b0\.\d+\b")
EFFECT = re.compile(r"surviv|viabilit|clonogen|cytotox|cell death|reduc|decreas|inhibit|kill", re.I)
# and, for treatment context, some marker that the NPI alone is meant
ALONE = re.compile(r"\balone\b|\bonly\b|untreated|control|monotherap|single agent|per se", re.I)


def get(url: str, timeout: int = 60) -> bytes:
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read()


def esearch(term: str, retmax: int) -> list[str]:
    q = urllib.parse.urlencode({"db": "pmc", "term": term, "retmax": retmax, "retmode": "json"})
    d = json.loads(get(f"{EUTILS}/esearch.fcgi?{q}"))
    return d.get("esearchresult", {}).get("idlist", [])


def efetch_text(ids: list[str]) -> dict[str, str]:
    """PMC id -> plain text of the article XML."""
    q = urllib.parse.urlencode({"db": "pmc", "id": ",".join(ids), "retmode": "xml"})
    raw = get(f"{EUTILS}/efetch.fcgi?{q}").decode("utf-8", "replace")
    out = {}
    for art in re.split(r"(?=<article[ >])", raw):
        m = re.search(r'pub-id-type="pmcid?"[^>]*>(?:PMC)?(\d+)<', art) or \
            re.search(r"<article-id[^>]*>PMC(\d+)</article-id>", art)
        if not m:
            continue
        txt = re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", art)))
        out[f"PMC{m.group(1)}"] = txt
    return out


def sentences(text: str) -> list[str]:
    return re.split(r"(?<=[.])\s+(?=[A-Z(])", text)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--npi", choices=sorted(QUERIES), help="use the built-in query set for this NPI")
    ap.add_argument("--query", action="append", help="extra PMC query (repeatable)")
    ap.add_argument("--cells", nargs="*", help="restrict to these curated cell-line ids or aliases")
    ap.add_argument("--max-papers", type=int, default=40, help="per query")
    ap.add_argument("--batch", type=int, default=15, help="papers per efetch call")
    ap.add_argument("--require-alone", action="store_true",
                    help="only sentences that also mention alone/control/untreated")
    ap.add_argument("--out", default="out/npi_monotherapy_candidates")
    a = ap.parse_args(argv)
    terms = list(a.query or [])
    if a.npi:
        terms += QUERIES[a.npi]
    if not terms:
        ap.error("pass --npi and/or --query")

    aliases = CELL_ALIASES
    if a.cells:
        want = {c.upper() for c in a.cells}
        aliases = {k: v for k, v in CELL_ALIASES.items()
                   if k.upper() in want or any(x.upper() in want for x in v)}
        if not aliases:
            aliases = {c: [c] for c in a.cells}
    pat = {cl: re.compile("|".join(re.escape(x) for x in names)) for cl, names in aliases.items()}

    ids: list[str] = []
    for t in terms:
        found = esearch(t, a.max_papers)
        print(f"  {len(found)} papers for: {t[:80]}...", file=sys.stderr)
        ids += found
        time.sleep(0.4)
    ids = list(dict.fromkeys(ids))
    print(f"{len(ids)} unique open-access papers", file=sys.stderr)

    rows = []
    for i in range(0, len(ids), a.batch):
        chunk = ids[i:i + a.batch]
        try:
            texts = efetch_text(chunk)
        except Exception as e:
            print(f"  efetch failed for {chunk[0]}...: {e}", file=sys.stderr)
            continue
        for pmcid, txt in texts.items():
            title = txt[:220].strip()
            for s in sentences(txt):
                if len(s) > 700 or not NUMBER.search(s) or not EFFECT.search(s):
                    continue
                if a.require_alone and not ALONE.search(s):
                    continue
                hits = [cl for cl, p in pat.items() if p.search(s)]
                if hits:
                    rows.append({"pmcid": pmcid, "cell_lines": ",".join(hits),
                                 "sentence": s.strip()[:600], "title_head": title[:120]})
        print(f"  scanned {min(i + a.batch, len(ids))}/{len(ids)}", file=sys.stderr)
        time.sleep(0.4)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    tag = a.npi or "custom"
    import pandas as pd

    t = pd.DataFrame(rows).drop_duplicates(subset=["pmcid", "sentence"])
    t.to_csv(out / f"{tag}.tsv", sep="\t", index=False)
    print(f"\n{len(t)} candidate sentences across {t['pmcid'].nunique() if len(t) else 0} papers "
          f"-> {out / f'{tag}.tsv'}")
    if len(t):
        print(t["cell_lines"].str.split(",").explode().value_counts().to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
