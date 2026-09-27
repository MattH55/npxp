"""Persistence: signature libraries as Parquet, patient states as .npz."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from .model import PatientState, Signature

_SCALARS = [
    "sig_id", "kind", "provenance", "modality", "tissue", "species", "duration",
    "intensity", "sample_size", "contrast", "quality_flag",
]


def save_signatures(sigs: list[Signature], path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    cols: dict[str, list] = {k: [getattr(s, k) for s in sigs] for k in _SCALARS}
    cols["source_accessions"] = [list(s.source_accessions) for s in sigs]
    cols["up_genes"] = [list(s.up_genes) for s in sigs]
    cols["down_genes"] = [list(s.down_genes) for s in sigs]
    cols["genes"] = [list(s.genes) for s in sigs]
    cols["z_vector"] = [s.z.astype(np.float32).tolist() for s in sigs]
    cols["meta_json"] = [json.dumps(s.meta, sort_keys=True, default=str) for s in sigs]
    schema = pa.schema(
        [(k, pa.string()) for k in _SCALARS if k != "sample_size"]
        + [("sample_size", pa.int64())]
        + [(k, pa.list_(pa.string())) for k in ("source_accessions", "up_genes", "down_genes", "genes")]
        + [("z_vector", pa.list_(pa.float32())), ("meta_json", pa.string())]
    )
    table = pa.table({f.name: cols[f.name] for f in schema}, schema=schema)
    pq.write_table(table, path)
    return path


def load_signatures(path: str | Path) -> list[Signature]:
    rows = pq.read_table(path).to_pylist()
    out = []
    for r in rows:
        out.append(
            Signature(
                **{k: r[k] for k in _SCALARS},
                genes=r["genes"],
                z=np.asarray(r["z_vector"], dtype=float),
                source_accessions=r["source_accessions"] or [],
                up_genes=r["up_genes"] or [],
                down_genes=r["down_genes"] or [],
                meta=json.loads(r["meta_json"] or "{}"),
            )
        )
    return out


def signature_index(sigs: list[Signature]) -> dict[str, Signature]:
    idx: dict[str, Signature] = {}
    for s in sigs:
        if s.sig_id in idx:
            raise ValueError(f"duplicate signature id {s.sig_id!r}")
        idx[s.sig_id] = s
    return idx


def save_patient(p: PatientState, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    meta = {
        "sample_id": p.sample_id, "tissue": p.tissue, "reference": p.reference,
        "flags": p.flags, "meta": p.meta,
    }
    with open(path, "wb") as fh:  # file handle stops numpy appending ".npz"
        np.savez_compressed(
            fh,
            genes=np.array(p.genes, dtype=str),
            expr_z=p.expr_z,
            disease_vector=p.disease_vector,
            meta_json=np.array(json.dumps(meta, default=str)),
        )
    return path


def load_patient(path: str | Path) -> PatientState:
    with np.load(path, allow_pickle=False) as d:
        meta = json.loads(str(d["meta_json"]))
        return PatientState(
            sample_id=meta["sample_id"],
            tissue=meta["tissue"],
            genes=[str(g) for g in d["genes"]],
            expr_z=d["expr_z"],
            disease_vector=d["disease_vector"],
            reference=meta["reference"],
            flags=meta.get("flags", []),
            meta=meta.get("meta", {}),
        )
