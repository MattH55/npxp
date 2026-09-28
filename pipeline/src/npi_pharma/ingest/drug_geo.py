"""Drug signatures from ordinary GEO series, for drugs LINCS does not cover.

LINCS Phase II holds only 5 of this project's 13 curated drugs and no platinum
agent, so cisplatin, oxaliplatin, carboplatin, 5-fluorouracil and erastin
signatures come from `configs/drug_signatures_geo.yaml`: acute drug-treated versus
vehicle experiments in human cancer cell lines, with the sample columns read off
GEO's own annotation.

Each record becomes a :class:`Signature` of kind DRUG, built by the same
unpaired moderated-DE path used for NPI signatures, so the two sides of a
retrieval comparison are processed identically.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from ..model import DRUG, Signature
from ..signatures.build import build_signature

PROV_GEO_DRUG = "geo_drug_treatment"


def _read_table(path: Path, fmt: dict[str, Any]) -> pd.DataFrame:
    """Read a processed expression table, skipping any leading comment lines."""
    enc = fmt.get("encoding", "utf-8")
    opener = (lambda: gzip.open(path, "rt", encoding=enc, errors="replace")) if path.suffix == ".gz" \
        else (lambda: open(path, encoding=enc, errors="replace"))
    skip = 0
    with opener() as fh:
        for line in fh:
            if line.startswith("#"):
                skip += 1
                continue
            break
    return pd.read_csv(path, sep=fmt.get("sep", "\t"), skiprows=skip, encoding=enc,
                       engine="python", on_bad_lines="skip")


def _gene_index(df: pd.DataFrame, fmt: dict[str, Any]) -> pd.Series:
    """Gene symbols for each row, from a symbol column or split out of a compound id."""
    if fmt.get("gene_col") and fmt["gene_col"] in df.columns:
        return df[fmt["gene_col"]].astype(str)
    sep = fmt.get("gene_from_id")
    idc = fmt.get("id_col")
    if sep and idc and idc in df.columns:
        # rows like ENSG00000000003_TSPAN6 -> TSPAN6
        return df[idc].astype(str).str.split(sep).str[-1]
    raise ValueError(f"cannot find gene symbols: format {fmt}")


def load_catalog(path: str | Path) -> list[dict[str, Any]]:
    return (yaml.safe_load(Path(path).read_text()) or {}).get("drug_signatures", [])


def build_from_record(rec: dict[str, Any], raw_dir: str | Path, min_genes: int = 2000) -> Signature:
    """Build one drug signature from a catalog record (control vs treated, unpaired)."""
    fmt = rec.get("format") or {}
    df = _read_table(Path(raw_dir) / rec["file"], fmt)
    df.columns = [str(c).strip() for c in df.columns]
    missing = [c for c in rec["control"] + rec["treated"] if c not in df.columns]
    if missing:
        raise ValueError(f"{rec['accession']} {rec['drug_id']}: columns not in file: {missing[:4]}")
    genes = _gene_index(df, fmt)
    expr = df[rec["control"] + rec["treated"]].apply(pd.to_numeric, errors="coerce")
    expr.index = genes.str.upper().str.strip()
    expr = expr[(expr.index != "") & ~expr.index.isin(["NAN", "NA", "NONE", "-"])]
    expr = expr.groupby(level=0).sum()                     # several rows per symbol -> sum
    # counts or FPKM both arrive as non-negative abundances; put on a log scale and
    # drop the low-abundance tail where log ratios are mostly noise
    total = expr.sum(axis=0).replace(0, np.nan)
    cpm = expr / total * 1e6
    keep = (cpm >= 1).mean(axis=1) >= 0.5
    logx = np.log2(cpm[keep] + 1)
    if len(logx) < min_genes:
        raise ValueError(f"{rec['accession']} {rec['drug_id']}: only {len(logx)} genes above the filter")
    sig_id = rec["drug_id"] if rec.get("cell_line") is None else f"{rec['drug_id']}|{rec['cell_line']}"
    meta = {
        "npi_id": sig_id,                                   # build_signature's id field
        "contrast": f"{rec.get('condition', 'drug')} vs control, unpaired",
        "provenance": PROV_GEO_DRUG,
        "source_accessions": [rec["accession"]],
        "notes": rec.get("notes"),
    }
    sig = build_signature(logx[rec["control"]], logx[rec["treated"]], meta, paired=False,
                          min_n=3)
    sig.kind = DRUG
    sig.tissue = rec.get("cell_line")
    sig.modality = "compound"
    sig.meta = dict(sig.meta or {}) | {
        "cell_line": rec.get("cell_line"), "cell_line_id": rec.get("cell_line_id"),
        "curated_cell_line": bool(rec.get("curated_cell_line")), "condition": rec.get("condition"),
        "drug_id": rec["drug_id"], "n_control": len(rec["control"]), "n_treated": len(rec["treated"]),
    }
    return sig


def build_all(catalog: str | Path, raw_dir: str | Path) -> tuple[list[Signature], list[str]]:
    sigs, problems = [], []
    for rec in load_catalog(catalog):
        try:
            sigs.append(build_from_record(rec, raw_dir))
        except (ValueError, KeyError, FileNotFoundError, OSError) as e:
            problems.append(f"{rec.get('accession')} {rec.get('drug_id')}: {e}")
    return sigs, problems
