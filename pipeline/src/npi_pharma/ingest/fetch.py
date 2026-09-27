"""Download GEO series matrices, platform annotations and LINCS files.

Everything comes from NCBI's public HTTPS mirror of the GEO FTP site
(https://ftp.ncbi.nlm.nih.gov/geo/), which needs no login. The host must be
reachable from wherever this runs.
"""

from __future__ import annotations

import gzip
import re
import shutil
import sys
import urllib.request
from pathlib import Path

import pandas as pd

GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo"


def _stem(acc: str) -> str:
    """GSE95640 -> GSE95nnn; GPL570 -> GPLnnn (GEO directory bucketing)."""
    prefix, digits = re.match(r"([A-Z]+)(\d+)$", acc).groups()
    return f"{prefix}{digits[:-3]}nnn" if len(digits) > 3 else f"{prefix}nnn"


def series_dir(gse: str, sub: str) -> str:
    return f"{GEO_FTP}/series/{_stem(gse)}/{gse}/{sub}/"


def platform_annot_url(gpl: str) -> str:
    return f"{GEO_FTP}/platforms/{_stem(gpl)}/{gpl}/annot/{gpl}.annot.gz"


def _get(url: str, timeout: int = 60) -> bytes:
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read()


def list_dir(url: str) -> list[str]:
    """File names linked from an NCBI FTP-over-HTTPS directory listing."""
    html = _get(url).decode("utf-8", "replace")
    return sorted({h for h in re.findall(r'href="([^"?/][^"]*)"', html) if not h.startswith("http")})


def download(url: str, dest: Path, overwrite: bool = False) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0 and not overwrite:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as fh:
        shutil.copyfileobj(r, fh, length=1 << 20)
    tmp.replace(dest)
    return dest


def series_platforms(matrix_path: Path) -> list[str]:
    with gzip.open(matrix_path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Series_platform_id"):
                return [p.strip('"') for p in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!series_matrix_table_begin"):
                break
    return []


def annot_to_probe_map(annot_gz: Path, out_tsv: Path) -> Path:
    """GPL .annot.gz (GEO's curated probe annotation) -> probe<TAB>symbol."""
    rows, header, inside = [], None, False
    with gzip.open(annot_gz, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                inside = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if inside:
                parts = line.rstrip("\n").split("\t")
                if header is None:
                    header = parts
                else:
                    rows.append(parts)
    if header is None or "Gene symbol" not in header:
        raise ValueError(f"{annot_gz}: no platform table with a 'Gene symbol' column")
    df = pd.DataFrame(rows, columns=header)[["ID", "Gene symbol"]]
    df = df[df["Gene symbol"].str.len() > 0]
    df.to_csv(out_tsv, sep="\t", header=False, index=False)
    return out_tsv


def fetch_geo_series(gse: str, raw_dir: str | Path, overwrite: bool = False, log=print) -> dict[str, list[Path]]:
    """Fetch all series matrices of a GSE plus a probe map per platform."""
    out = Path(raw_dir) / gse
    names = [n for n in list_dir(series_dir(gse, "matrix")) if n.endswith("_series_matrix.txt.gz")]
    if not names:
        raise FileNotFoundError(f"{gse}: no series matrix files listed")
    got: dict[str, list[Path]] = {"series_matrix": [], "probe_map": []}
    for n in names:
        p = download(series_dir(gse, "matrix") + n, out / n, overwrite)
        got["series_matrix"].append(p)
        log(f"  {p}")
        for gpl in series_platforms(p):
            pm = out / (f"probe_map_{gpl}.tsv" if len(names) > 1 else "probe_map.tsv")
            if pm.exists() and not overwrite:
                got["probe_map"].append(pm)
                continue
            try:
                annot = download(platform_annot_url(gpl), out / f"{gpl}.annot.gz", overwrite)
                got["probe_map"].append(annot_to_probe_map(annot, pm))
                log(f"  {pm} (from {gpl}.annot.gz)")
            except Exception as e:  # sequencing platforms have no .annot; matrices may already be symbols
                print(f"  no probe map for {gpl}: {e}", file=sys.stderr)
    return got


LINCS_PATTERNS = {
    "gctx": r"_Level5_COMPZ.*\.gctx\.gz$",
    "sig_info": r"_sig_info.*\.txt\.gz$",
    "gene_info": r"_gene_info.*\.txt\.gz$",
    "sig_metrics": r"_sig_metrics.*\.txt\.gz$",
    "pert_info": r"_pert_info.*\.txt\.gz$",
}


def fetch_lincs(gse: str, raw_dir: str | Path, overwrite: bool = False, log=print) -> dict[str, Path]:
    """Fetch and gunzip the LINCS Level 5 GCTX and its metadata from a GEO series
    (GSE70138 Phase II, ~2 GB compressed; GSE92742 Phase I is ~20 GB)."""
    out = Path(raw_dir) / gse
    names = list_dir(series_dir(gse, "suppl"))
    got = {}
    for key, pat in LINCS_PATTERNS.items():
        hits = [n for n in names if re.search(pat, n)]
        if not hits:
            if key in ("gctx", "sig_info", "gene_info"):
                raise FileNotFoundError(f"{gse}: no supplementary file matching {pat}")
            continue
        gz = download(series_dir(gse, "suppl") + hits[-1], out / hits[-1], overwrite)
        plain = gz.with_suffix("")
        if overwrite or not plain.exists():
            with gzip.open(gz, "rb") as src, open(plain, "wb") as dst:
                shutil.copyfileobj(src, dst, length=1 << 20)
        if key == "gctx":
            gz.unlink()  # keep only the uncompressed HDF5 to save disk
        got[key] = plain
        log(f"  {key}: {plain}")
    return got
