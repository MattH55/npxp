"""Rank NPI x drug pairs for one patient."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from ..config import canonical_tissue
from ..genes import GeneOverlapError
from ..model import PatientState, Signature
from .score import COMPONENTS, score_pair

TABLE_COLS = [
    "composite_rank", "npi_id", "drug_id", "composite", "composite_raw", "confidence",
    *COMPONENTS, "complementarity_gain", "orthogonality_raw", "reverse_npi", "reverse_drug", "reverse_combo", "n_genes", "flags",
]


def select_npis(npis: list[Signature], patient: PatientState, cfg: dict, npi_class: str | None = None,
                npi_ids: list[str] | None = None, allow_mismatch: bool = False) -> tuple[list[Signature], list[str]]:
    """Filter NPIs by class/id and tissue compatibility. Returns (kept, skip reasons)."""
    kept, skipped = [], []
    pt = canonical_tissue(patient.tissue, cfg)
    for s in npis:
        if npi_ids and s.sig_id not in npi_ids:
            continue
        if npi_class and (s.modality or "").lower() != npi_class.lower():
            continue
        if not allow_mismatch and not npi_ids and canonical_tissue(s.tissue, cfg) != pt:
            skipped.append(f"{s.sig_id}: tissue {s.tissue} incompatible with patient {patient.tissue}")
            continue
        kept.append(s)
    return kept, skipped


def rank_pairs(
    patient: PatientState,
    npis: list[Signature],
    drugs: list[Signature],
    gene_sets: dict[str, list[str]],
    cfg: dict[str, Any],
) -> tuple[pd.DataFrame, dict[tuple[str, str], dict[str, Any]]]:
    """Score every NPI x drug pair. Pairs that cannot be scored stay in the table
    with NaN scores and an error flag instead of disappearing."""
    rows, reports = [], {}
    for npi in npis:
        for drug in drugs:
            try:
                rep = score_pair(patient, npi, drug, gene_sets, cfg)
            except GeneOverlapError as e:
                rows.append({"npi_id": npi.sig_id, "drug_id": drug.sig_id, "flags": f"error:{e}"})
                continue
            reports[(npi.sig_id, drug.sig_id)] = rep
            rows.append({k: rep[k] for k in TABLE_COLS if k in rep} | {"flags": ";".join(rep["flags"])})
    df = pd.DataFrame(rows, columns=[c for c in TABLE_COLS if c != "composite_rank"])
    df = df.sort_values(["composite", "npi_id", "drug_id"], ascending=[False, True, True], na_position="last")
    df.insert(0, "composite_rank", np.where(df["composite"].notna(), np.arange(1, len(df) + 1), np.nan))
    df = df.reset_index(drop=True)
    for _, r in df.iterrows():
        rep = reports.get((r["npi_id"], r["drug_id"]))
        if rep is not None:
            rep["composite_rank"] = int(r["composite_rank"])
            rep["n_pairs_ranked"] = int(df["composite"].notna().sum())
    return df, reports
