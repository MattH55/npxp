"""Fetch the real external data the npi_pharma.mechanism module needs.

    python scripts/fetch_mechanism_data.py

Both sources are real, open, and ungated:
  - SynLethDB 3.0 (Feng, Zhang & Zheng, ShanghaiTech; CC-BY-4.0), the curated
    (non-computational) human synthetic-lethality table, from its Zenodo
    record (DOI 10.5281/zenodo.22843223) -- a stable, versioned archive
    rather than scraping the live site.
  - DGIdb 5.0's current interactions.tsv bulk export (dgidb.org/downloads).
"""
from __future__ import annotations

import os
import sys
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SYNLETHDB_URL = "https://zenodo.org/api/records/22843223/files/Human.SL.detailed.tsv/content"
SYNLETHDB_OUT = os.path.join(ROOT, "data", "raw", "synlethdb", "Human.SL.detailed.tsv")
DGIDB_URL = "https://dgidb.org/data/latest/interactions.tsv"
DGIDB_OUT = os.path.join(ROOT, "data", "raw", "dgidb", "interactions.tsv")


def fetch(url: str, dest: str) -> None:
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if os.path.isfile(dest):
        print(f"  already have {dest} ({os.path.getsize(dest):,} bytes)")
        return
    print(f"  downloading {url}")
    urllib.request.urlretrieve(url, dest)
    print(f"  wrote {dest} ({os.path.getsize(dest):,} bytes)")


def main() -> int:
    print("SynLethDB 3.0 (real, curated human SL pairs, CC-BY-4.0):")
    fetch(SYNLETHDB_URL, SYNLETHDB_OUT)
    print("DGIdb 5.0 (real drug-gene interaction claims):")
    fetch(DGIDB_URL, DGIDB_OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
