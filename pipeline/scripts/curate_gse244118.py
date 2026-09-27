"""Curate the GSE244118 sample sheet (metabolic group per counts column).

GSE244118 holds abdominal subcutaneous adipose RNA-seq from metabolically
healthy lean (MHL), metabolically healthy obese (MHO) and metabolically
unhealthy obese (MUO) adults. It serves as the lean-vs-obese contrast that
places an obese patient relative to health (`npi-pharma build-offset`).

The counts columns are library names, and GEO keeps the group in
`metabolic group`. The series-matrix title carries the same library name in
brackets. This script joins the two, checks that the group tag inside each
library name agrees with GEO's `metabolic group`, and writes
`samples_curated.tsv`.

Run after `npi-pharma fetch-geo GSE244118 --suppl`:

    python scripts/curate_gse244118.py --raw-dir data/raw
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from npi_pharma.ingest.geo import read_counts, read_series_matrix

GSE = "GSE244118"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--raw-dir", default="data/raw")
    a = ap.parse_args(argv)
    d = Path(a.raw_dir) / GSE

    _, sm = read_series_matrix(d / f"{GSE}_series_matrix.txt.gz")
    sm["sample"] = sm["title"].str.extract(r"\[(.+)\]", expand=False)
    cols = set(read_counts(d / f"{GSE}_abdominal.fat_all.gene_counts.txt.gz").columns)
    missing = sorted(set(sm["sample"]) - cols)
    tag = sm["sample"].str.extract(r"_(MHL|MHO|MUO)_", expand=False)
    disagree = int((tag != sm["metabolic group"]).sum())
    print(f"samples: {len(sm)}; not in counts: {len(missing)}; group tag disagrees with GEO: {disagree}")
    print(sm["metabolic group"].value_counts().to_string())
    if missing or disagree:
        print("checks failed; not writing the curated sheet", file=sys.stderr)
        return 1
    # Library-name prefix = sequencing project. All 4 "PE" libraries are MUO, so
    # keep it as `batch` and build the offset from one batch (--filter batch=PSQ).
    sm["batch"] = sm["sample"].str.extract(r"^([A-Za-z]+)", expand=False)
    print(sm.groupby(["metabolic group", "batch"]).size().to_string())
    out = d / "samples_curated.tsv"
    sm.rename(columns={"metabolic group": "group"})[["sample", "group", "batch", "geo_accession"]].to_csv(
        out, sep="\t", index=False)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
