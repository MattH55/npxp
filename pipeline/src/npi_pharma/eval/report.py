"""JSON / TSV reports."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .. import EVIDENCE_TIER, __version__

DISCLAIMER = (
    "Tier 3 (mechanism-only) signature-composition scores. These rank how "
    "plausibly an NPI and a drug complement each other against this "
    "transcriptome; they are not measured or calibrated synergy, not a PK "
    "prediction, and not clinical advice."
)


def _default(o: Any):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))


def explain(rep: dict[str, Any], n_pathways: int = 5, n_genes: int = 10) -> dict[str, Any]:
    """Short explanation of what drove a pair's score."""
    pw = rep["pathways"][:n_pathways]
    return {
        "pair": f"{rep['npi_id']} x {rep['drug_id']}",
        "composite": rep["composite"],
        "composite_rank": rep.get("composite_rank"),
        "summary": (
            f"combo reversal {rep['reverse_combo']:.3f} vs best single agent "
            f"{max(rep['reverse_npi'], rep['reverse_drug']):.3f} (gain {rep['complementarity_gain']:+.3f}, {rep['complementarity']:+.2f} of headroom); "
            f"pathway_joint {rep['pathway_joint']:+.3f}; confidence x{rep['confidence']:.2f}"
        ),
        "pathways_driving": [
            {k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()} for r in pw
        ],
        "top_genes": [g["gene"] for g in rep["genes_driving_reversal"][:n_genes]],
        "flags": rep["flags"],
        "warnings": rep["warnings"],
    }


def envelope(body: dict[str, Any], inputs: dict[str, Any]) -> dict[str, Any]:
    return {
        "tool": "npi-pharma", "version": __version__, "evidence_tier": EVIDENCE_TIER,
        "disclaimer": DISCLAIMER, "inputs": inputs, **body,
    }


def _sanitize(o: Any) -> Any:
    if isinstance(o, dict):
        return {str(k): _sanitize(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_sanitize(v) for v in o]
    if isinstance(o, (float, np.floating)):
        return None if np.isnan(o) else float(o)
    return o


def write_json(obj: dict[str, Any], path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_sanitize(obj), indent=2, default=_default, allow_nan=False) + "\n")
    return path


def write_tsv(df: pd.DataFrame, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, float_format="%.6f")
    return path
