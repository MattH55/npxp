"""What LINCS Phase I (GSE92742) would add, from its metadata alone.

Phase II (GSE70138) lacks the platinums, 5-FU and erastin, which is why this project
built GEO consensus signatures instead -- and that route failed: no drug reached the
0.30 cross-series agreement the validation set for "medium", and no platinum reached
even 0.10 (docs/drug_consensus.md). Phase I is the obvious next place to look.

Its Level 5 matrix is 20.3 GB **gzipped**, so before spending that this asks the
decision-relevant question from the 12 MB of metadata: which drugs does Phase I
actually carry, and in how many cell lines? Cell-line count is the thing that
matters, because a LINCS consensus is a median over cell lines and that is precisely
what makes it trustworthy -- the same drug in two cell lines agrees at only 0.18
(docs/npi_drug_retrieval.md).

It also reports whether the matrix could be fetched here at all. A .gctx is HDF5,
which needs random access, and a gzip stream cannot be seeked -- so there is no way
to read a few columns out of the remote file. It is all 20.3 GB or nothing.

    python scripts/scope_lincs_phase1.py
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import urllib.request
from pathlib import Path

import pandas as pd

BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/"
METADATA = ["GSE92742_Broad_LINCS_sig_info.txt.gz",
            "GSE92742_Broad_LINCS_pert_info.txt.gz",
            "GSE92742_Broad_LINCS_gene_info.txt.gz"]
LEVEL5 = "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz"

# The drugs this project needs and Phase II does not carry, plus the two it does
# (fluorouracil, mitomycin-c) as a check that the lookup finds what is really there.
WANTED = ["cisplatin", "carboplatin", "oxaliplatin", "fluorouracil", "erastin",
          "triapine", "lomustine", "cyclophosphamide", "methotrexate", "pemetrexed",
          "topotecan", "idarubicin", "vincristine", "vinblastine", "vinorelbine",
          "melphalan", "carmustine", "bleomycin", "actinomycin", "azacitidine",
          "mitomycin", "dactinomycin"]

# Cell-line counts that would earn each reliability band. A LINCS consensus is "high"
# in this project because Phase II medians over 7 cell lines per compound.
BANDS = [(7, "high -- at or above the Phase II median of 7 cell lines"),
         (4, "usable -- fewer cell lines than Phase II typically has"),
         (2, "weak -- a median over 2 or 3 lines barely averages anything"),
         (1, "useless -- a single cell line is not a consensus")]


def band_for(n_cell_lines: int) -> str:
    for threshold, label in BANDS:
        if n_cell_lines >= threshold:
            return label
    return "absent"


def remote_size_mb(name: str) -> float | None:
    try:
        req = urllib.request.Request(BASE + name, method="HEAD")
        with urllib.request.urlopen(req, timeout=90) as h:
            return int(h.headers["Content-Length"]) / 2**20
    except Exception as e:
        print(f"  (could not size {name}: {e})", file=sys.stderr)
        return None


def fetch_metadata(out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    for n in METADATA:
        d = out / n
        if d.exists() and d.stat().st_size > 0:
            continue
        with urllib.request.urlopen(BASE + n, timeout=300) as r, open(d, "wb") as fh:
            shutil.copyfileobj(r, fh, length=1 << 20)
        print(f"  fetched {n} ({d.stat().st_size / 2**20:.1f} MB)", file=sys.stderr)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--raw-dir", default="data/raw/GSE92742")
    ap.add_argument("--phase2", default="data/processed/signatures/lincs_all.parquet")
    ap.add_argument("--out", default="out/lincs_phase1_scope")
    a = ap.parse_args(argv)
    raw, out = Path(a.raw_dir), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    fetch_metadata(raw)

    si = pd.read_csv(raw / METADATA[0], sep="\t", low_memory=False)
    si = si[si["pert_type"] == "trt_cp"]
    name = si["pert_iname"].astype(str).str.lower()

    have2: set[str] = set()
    if Path(a.phase2).exists():
        sys.path.insert(0, str(Path(__file__).parents[1] / "src"))
        from npi_pharma.store import load_signatures
        have2 = {s.sig_id.lower() for s in load_signatures(a.phase2)}

    rows = []
    for w in WANTED:
        hit = si[name.str.contains(w, regex=False, na=False)]
        if not len(hit):
            rows.append({"query": w, "pert_iname": None, "n_signatures": 0,
                         "n_cell_lines": 0, "band": "absent", "in_phase2": w in have2})
            continue
        for pn, g in hit.groupby("pert_iname"):
            rows.append({"query": w, "pert_iname": pn, "n_signatures": len(g),
                         "n_cell_lines": int(g["cell_id"].nunique()),
                         "n_doses": int(g["pert_idose"].nunique()),
                         "band": band_for(int(g["cell_id"].nunique())),
                         "in_phase2": str(pn).lower() in have2})
    t = pd.DataFrame(rows).sort_values(["n_cell_lines", "n_signatures"], ascending=False)
    t.to_csv(out / "coverage.tsv", sep="\t", index=False)

    print(f"\nLINCS Phase I coverage for the {len(WANTED)} drugs this project needs\n")
    print(t[["query", "pert_iname", "n_signatures", "n_cell_lines", "band",
             "in_phase2"]].to_string(index=False, na_rep="-"))

    gained = t[(~t["in_phase2"]) & (t["n_cell_lines"] >= 7)]
    absent = t[t["n_cell_lines"] == 0]
    size = remote_size_mb(LEVEL5)
    free = shutil.disk_usage(".").free / 2**20
    verdict = {
        "would_gain_high_reliability": sorted(gained["pert_iname"].dropna().unique()),
        "absent_from_phase1": sorted(absent["query"].unique()),
        "level5_gz_mb": round(size, 1) if size else None,
        "free_mb": round(free, 1),
        "fits": bool(size and free > size),
        "why_no_partial_fetch":
            "A .gctx is HDF5 and needs random access; a gzip stream cannot be seeked, "
            "so no subset of the remote matrix can be read. It is all of it or none.",
        "decompressed_also_needed":
            "Reading it requires the decompressed HDF5 as well, which is several times "
            "the 20.3 GB archive, so the archive size is a floor on the space needed, "
            "not the total.",
    }
    (out / "verdict.json").write_text(json.dumps(verdict, indent=2))

    print(f"\nwould reach 'high' reliability and are missing today: "
          f"{', '.join(verdict['would_gain_high_reliability']) or 'none'}")
    print(f"absent from Phase I entirely: {', '.join(verdict['absent_from_phase1']) or 'none'}")
    if size:
        print(f"\nLevel 5 matrix: {size / 1024:.1f} GB gzipped; {free / 1024:.1f} GB free here "
              f"-> {'fits' if verdict['fits'] else 'DOES NOT FIT'}")
        print(verdict["why_no_partial_fetch"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
