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


class TooLarge(Exception):
    """A download exceeded its size cap and was abandoned."""


def download(url: str, dest: Path, overwrite: bool = False,
             max_bytes: int | None = None) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0 and not overwrite:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=120) as r:
        if max_bytes is not None:
            declared = r.headers.get("Content-Length")
            if declared is not None and int(declared) > max_bytes:
                raise TooLarge(f"{url}: {int(declared) / 2**20:.0f} MB "
                               f"exceeds the {max_bytes / 2**20:.0f} MB cap")
        with open(tmp, "wb") as fh:
            if max_bytes is None:
                shutil.copyfileobj(r, fh, length=1 << 20)
            else:
                # Some listings declare no length, so cap the stream as well and drop
                # the partial file: a single unbounded supplementary file (a 2.7 GB
                # Hi-C map shipped beside a count table) can fill the disk and take
                # the whole run down with it.
                written = 0
                while chunk := r.read(1 << 20):
                    written += len(chunk)
                    if written > max_bytes:
                        fh.close()
                        tmp.unlink(missing_ok=True)
                        raise TooLarge(f"{url}: exceeds the {max_bytes / 2**20:.0f} MB cap")
                    fh.write(chunk)
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


PLATFORM_TABLE_URL = ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
                      "?acc={gpl}&targ=self&view=data&form=text")
SYMBOL_COL = re.compile(r"^(gene[ _]?symbol|symbol|gene[ _]?name)$", re.I)


def probe_map_from_platform_table(gpl: str, out_tsv: Path, timeout: int = 180) -> Path:
    """probe<TAB>symbol from GEO's platform table, for platforms with no ``.annot.gz``.

    Many newer or vendor-specific platforms (Affymetrix PrimeView, Agilent arrays)
    have no curated ``.annot`` file, but their submitted platform table does carry a
    gene-symbol column. Raises if the table has no such column (some platforms give
    only GenBank accessions).
    """
    raw = _get(PLATFORM_TABLE_URL.format(gpl=gpl), timeout).decode("utf-8", "replace")
    header, rows = None, []
    for line in raw.splitlines():
        if line.startswith(("^", "!", "#")) or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if header is None:
            if parts[0].strip().upper() != "ID":
                continue
            header = [p.strip() for p in parts]
            continue
        rows.append(parts)
    if header is None:
        raise ValueError(f"{gpl}: no platform table returned")
    sym = next((i for i, h in enumerate(header) if SYMBOL_COL.match(h)), None)
    assign = next((i for i, h in enumerate(header) if h.strip().lower() == "gene_assignment"), None)
    if sym is None and assign is None:
        raise ValueError(f"{gpl}: platform table has no gene-symbol column (has {header[:8]})")
    out = []
    for r in rows:
        if not r or not r[0].strip():
            continue
        if sym is not None and len(r) > sym and r[sym].strip():
            out.append((r[0].strip(), r[sym].strip().split("///")[0].strip()))
        elif assign is not None and len(r) > assign:
            # Affymetrix gene_assignment: "NM_001005484 // SAMD11 // description // ..."
            parts = [p.strip() for p in r[assign].split("///")[0].split("//")]
            if len(parts) >= 2 and parts[1]:
                out.append((r[0].strip(), parts[1]))
    if not out:
        raise ValueError(f"{gpl}: platform table gave no probe-to-symbol rows")
    pd.DataFrame(out).to_csv(out_tsv, sep="\t", header=False, index=False)
    return out_tsv


# Supplementary files this package can actually read as a table: a delimited text
# file or a spreadsheet, optionally compressed. Everything else a series ships
# beside its counts -- Hi-C contact maps (.hic), 10x Loupe projects (.cloupe),
# alignments, coverage tracks, raw-read tarballs, images -- cannot be parsed by
# supplementary_table() and runs to gigabytes, so it is never worth fetching.
SUPPL_TABLE = re.compile(r"\.(txt|tsv|csv|tab|xlsx?|mtx)(\.(gz|bz2|xz|zip))?$", re.I)
MAX_SUPPL_MB = 512


def is_suppl_table(name: str) -> bool:
    """Could this supplementary file name hold a readable expression table?"""
    return bool(SUPPL_TABLE.search(name)) and ".tar" not in name.lower()


def fetch_geo_series(
    gse: str, raw_dir: str | Path, overwrite: bool = False, suppl: bool = False, log=print,
    max_suppl_mb: int = MAX_SUPPL_MB,
) -> dict[str, list[Path]]:
    """Fetch all series matrices of a GSE plus a probe map per platform.

    With ``suppl`` also fetch the series' supplementary files. RNA-seq series keep
    their counts there; the series matrix has no table. Only files that could hold a
    readable table are fetched (:func:`is_suppl_table`), and each is capped at
    ``max_suppl_mb``; the rest are listed under ``"suppl_skipped"``. Without those
    two limits one series can exhaust the disk -- GSE236253 ships a 2.7 GB Hi-C map
    next to its 4.8 MB count table, and GSE276609 ships 8.5 GB of Loupe projects.
    """
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
            except Exception as annot_err:
                # No curated .annot: fall back to the submitted platform table, which
                # often carries a gene-symbol column. Sequencing platforms have neither,
                # and their matrices are usually already keyed by symbol or Ensembl id.
                try:
                    got["probe_map"].append(probe_map_from_platform_table(gpl, pm))
                    log(f"  {pm} (from the {gpl} platform table)")
                except Exception as table_err:
                    print(f"  no probe map for {gpl}: {annot_err}; platform table: {table_err}",
                          file=sys.stderr)
    if suppl:
        got["suppl"] = []
        got["suppl_skipped"] = []
        for n in list_dir(series_dir(gse, "suppl")):
            if not n.startswith(gse) or n.endswith("_RAW.tar"):
                continue
            if not is_suppl_table(n):
                got["suppl_skipped"].append(n)
                continue
            try:
                got["suppl"].append(download(series_dir(gse, "suppl") + n, out / n, overwrite,
                                             max_bytes=max_suppl_mb * 2**20))
                log(f"  {got['suppl'][-1]}")
            except TooLarge as e:
                got["suppl_skipped"].append(n)
                print(f"  skipped {n}: {e}", file=sys.stderr)
    return got


GENE_INFO_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"


def gene_info_to_id_map(gene_info_gz: Path, out_tsv: Path) -> Path:
    """NCBI gene_info -> ``ensembl_id<TAB>symbol`` (unversioned ENSG ids).

    An Ensembl id cross-referenced by several NCBI genes keeps the first
    protein-coding one, else the first listed.
    """
    gi = pd.read_csv(gene_info_gz, sep="\t", usecols=["Symbol", "dbXrefs", "type_of_gene"], dtype=str)
    gi["ens"] = gi["dbXrefs"].str.findall(r"Ensembl:(ENSG\d+)")
    gi = gi.explode("ens").dropna(subset=["ens"])
    gi["_pc"] = gi["type_of_gene"] != "protein-coding"
    gi = gi.sort_values("_pc", kind="stable").drop_duplicates("ens")
    gi[["ens", "Symbol"]].sort_values("ens").to_csv(out_tsv, sep="\t", header=False, index=False)
    return out_tsv


def fetch_gene_info(raw_dir: str | Path, overwrite: bool = False, log=print) -> Path:
    """Download NCBI Homo_sapiens.gene_info and write ``ensembl_to_symbol.tsv``."""
    raw_dir = Path(raw_dir)
    gz = download(GENE_INFO_URL, raw_dir / "Homo_sapiens.gene_info.gz", overwrite)
    out = gene_info_to_id_map(gz, raw_dir / "ensembl_to_symbol.tsv")
    log(f"  {out}")
    return out


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
