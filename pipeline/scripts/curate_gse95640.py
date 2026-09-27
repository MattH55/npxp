"""Curate the GSE95640 sample sheet (subject + timepoint) for `ingest-npi`.

GEO gives no subject id for GSE95640 (382 adipose RNA-seq samples, 191 obese
non-diabetic subjects, baseline vs after an 8-week 800-1000 kcal/d LCD). The
supplementary sample sheet lists the samples in subject pairs. This script turns
that adjacency into subject ids and refuses to write the sheet unless the
pairing and the timepoint direction both check out:

1. every adjacent pair has one CID1 and one CID2 sample, and the same sex and
   age in the series-matrix characteristics;
2. within-pair expression correlation (gene-centred log-CPM) clearly exceeds
   the between-subject correlation;
3. CID1 is the baseline. The series matrix says "CID1 is for the baseline
   without LCD / CID2 is after 8 weeks of LCD". The supplementary sheet's
   header comment says "CID1: after 8-week LCD / CID2: 6-month after LCD",
   which contradicts it. The data side with the series matrix: SCD, the
   canonical LCD-repressed adipose lipogenesis gene, must fall from CID1 to CID2.

Run after `npi-pharma fetch-geo GSE95640 --suppl --gene-info`:

    python scripts/curate_gse95640.py --raw-dir data/raw
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from npi_pharma.ingest.geo import counts_to_log_cpm, read_counts, read_series_matrix

GSE = "GSE95640"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--raw-dir", default="data/raw")
    a = ap.parse_args(argv)
    d = Path(a.raw_dir) / GSE

    sheet = pd.read_csv(d / f"{GSE}_sample_cid12.tsv.gz", sep="\t", comment="#", dtype=str)
    _, sm = read_series_matrix(d / f"{GSE}_series_matrix.txt.gz")
    sm = sm.set_index("title")
    ids = sheet["RNASEQ_SAMPLEID"].tolist()
    if len(ids) % 2 or set(ids) - set(sm.index):
        raise SystemExit("sample sheet does not match the series matrix titles")

    sheet["subject"] = [f"S{i // 2 + 1:03d}" for i in range(len(ids))]
    sheet = sheet.join(sm[["geo_accession", "gender", "age", "time"]], on="RNASEQ_SAMPLEID")
    if not (sheet["time"].str[:4] == sheet["TIMEPOINT"]).all():
        raise SystemExit("series-matrix time and sample-sheet TIMEPOINT disagree")
    g = sheet.groupby("subject")
    bad = [s for s, grp in g if sorted(grp["TIMEPOINT"]) != ["CID1", "CID2"]
           or grp["gender"].nunique() > 1 or grp["age"].nunique() > 1]
    print(f"pairs: {g.ngroups}; inconsistent sex/age/timepoint: {len(bad)}")

    idm = pd.read_csv(Path(a.raw_dir) / "ensembl_to_symbol.tsv", sep="\t", header=None, index_col=0, dtype=str)[1]
    expr = counts_to_log_cpm(read_counts(d / f"{GSE}_cid12_cng_rawcnts_qced.tsv.gz"), idm)
    z = expr[ids].sub(expr[ids].mean(axis=1), axis=0).to_numpy()
    corr = np.corrcoef(z.T)
    a_idx, b_idx = np.arange(0, len(ids), 2), np.arange(1, len(ids), 2)
    within = corr[a_idx, b_idx]
    between = corr[a_idx, np.roll(b_idx, 1)]
    print(f"within-pair r median {np.median(within):.3f}; shifted-pair r median {np.median(between):.3f}")

    t = sheet.set_index("RNASEQ_SAMPLEID")
    pre = t.index[t["TIMEPOINT"] == "CID1"]
    post = t.index[t["TIMEPOINT"] == "CID2"]
    pre = expr[sorted(pre, key=lambda s: t.loc[s, "subject"])]
    post = expr[sorted(post, key=lambda s: t.loc[s, "subject"])]
    scd = float((post.loc["SCD"].to_numpy() - pre.loc["SCD"].to_numpy()).mean())
    print(f"SCD mean log2FC CID2 - CID1: {scd:+.2f}")

    if bad or np.median(within) < np.median(between) + 0.2 or scd > -0.5:
        print("checks failed; not writing the curated sheet", file=sys.stderr)
        return 1
    out = d / "samples_curated.tsv"
    sheet.rename(columns={"RNASEQ_SAMPLEID": "sample", "TIMEPOINT": "timepoint", "gender": "sex"})[
        ["sample", "subject", "timepoint", "sex", "age", "geo_accession"]
    ].to_csv(out, sep="\t", index=False)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
