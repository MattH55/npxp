"""NPI catalog: one record per study / arm / tissue, validated on load."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from ..signatures.build import build_signature
from ..model import PROV_GEO, PROV_LOCAL, Signature
from . import geo

REQUIRED = ("npi_id", "modality", "tissue", "species", "contrast", "source_accessions")
METADATA = ("duration", "intensity", "sample_size")


@dataclass
class CatalogEntry:
    record: dict[str, Any]
    issues: list[str] = field(default_factory=list)

    @property
    def npi_id(self) -> str:
        return self.record["npi_id"]

    @property
    def quality_flag(self) -> str:
        r = self.record
        if r.get("exploratory"):
            return "exploratory"
        n = r.get("sample_size")
        if n is not None and n < r.get("_min_n", 6):
            return "small_n"
        if any(r.get(k) is None for k in METADATA):
            return "unverified_metadata"
        return "ok"


def load_catalog(path: str | Path, min_n: int = 6) -> list[CatalogEntry]:
    data = yaml.safe_load(Path(path).read_text()) or {}
    entries, seen = [], set()
    for rec in data.get("npis", []):
        issues = [f"missing required field {k!r}" for k in REQUIRED if not rec.get(k)]
        if rec.get("npi_id") in seen:
            raise ValueError(f"duplicate npi_id {rec['npi_id']!r} in {path}")
        seen.add(rec.get("npi_id"))
        issues += [f"{k} not yet curated" for k in METADATA if rec.get(k) is None]
        rec["_min_n"] = min_n
        entries.append(CatalogEntry(rec, issues))
    return entries


def _metadata(rec: dict[str, Any], provenance: str, flag: str) -> dict[str, Any]:
    md = {k: v for k, v in rec.items() if not k.startswith("_")}
    md["provenance"] = provenance
    md["quality_flag"] = flag
    return md


def build_from_entry(entry: CatalogEntry, raw_dir: str | Path) -> Signature:
    """Build a signature from the entry's ``inputs`` block.

    Supported inputs (paths relative to ``raw_dir``):
      series_matrix + probe_map (TSV: probe<TAB>symbol) + pairing
      expression (genes x samples TSV) + samples (TSV, first column = sample id) + pairing
      counts (raw RNA-seq counts, CSV/TSV) + id_map (TSV: gene id<TAB>symbol)
        + samples + pairing; converted to log2(CPM+1) with a low-expression filter
    """
    rec, raw_dir = entry.record, Path(raw_dir)
    inp = rec.get("inputs") or {}
    pairing = inp.get("pairing") or {}
    if not pairing.get("time_field") or pairing.get("pre") is None or pairing.get("post") is None:
        raise ValueError(f"{entry.npi_id}: inputs.pairing (time_field/pre/post) not curated")
    if "series_matrix" in inp:
        expr, samples = geo.read_series_matrix(raw_dir / inp["series_matrix"])
        if inp.get("probe_map"):
            pm = pd.read_csv(raw_dir / inp["probe_map"], sep="\t", header=None, index_col=0, dtype=str)[1]
            expr = geo.collapse_probes(expr, pm)
        prov = PROV_GEO
    elif "counts" in inp:
        idm = None
        if inp.get("id_map"):
            idm = pd.read_csv(raw_dir / inp["id_map"], sep="\t", header=None, index_col=0, dtype=str)[1]
        expr = geo.counts_to_log_cpm(geo.read_counts(raw_dir / inp["counts"]), idm)
        samples = pd.read_csv(raw_dir / inp["samples"], sep="\t", index_col=0, dtype=str)
        prov = PROV_GEO
    elif "expression" in inp:
        expr = pd.read_csv(raw_dir / inp["expression"], sep="\t", index_col=0)
        samples = pd.read_csv(raw_dir / inp["samples"], sep="\t", index_col=0, dtype=str)
        prov = inp.get("provenance", PROV_LOCAL)
    else:
        raise ValueError(f"{entry.npi_id}: no inputs.series_matrix, inputs.counts or inputs.expression")
    expr = geo.maybe_log2(expr)
    pre, post, paired = geo.split_pre_post(expr, samples, pairing)
    md = _metadata(rec, prov, entry.quality_flag)
    md["pairing"] = pairing
    if md.get("sample_size") is None:
        md["sample_size"] = int(min(pre.shape[1], post.shape[1]))
        if md["sample_size"] < rec["_min_n"] and md["quality_flag"] != "exploratory":
            md["quality_flag"] = "small_n"
    return build_signature(pre, post, md, paired=paired, min_n=rec["_min_n"])
