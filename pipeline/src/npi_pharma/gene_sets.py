"""Gene sets and pathway-activity scoring (shared by interventions and patients)."""

from __future__ import annotations

from importlib import resources
from pathlib import Path

import numpy as np

from .genes import normalize_symbol

DEFAULT_GMT = "curated_core.gmt"


def load_gmt(path: str | Path | None = None) -> dict[str, list[str]]:
    """Load a GMT file. ``None`` loads the bundled hand-curated core sets.

    The bundled sets are compact, hand-curated cores (not MSigDB). For real
    analyses pass the MSigDB Hallmark/KEGG GMT via ``--gene-sets``.
    """
    if path is None:
        text = resources.files("npi_pharma").joinpath("data", DEFAULT_GMT).read_text()
    else:
        text = Path(path).read_text()
    sets: dict[str, list[str]] = {}
    for line in text.splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        sets[parts[0]] = list(dict.fromkeys(normalize_symbol(g) for g in parts[2:] if g))
    return sets


def pathway_activity(
    genes: list[str], z: np.ndarray, gene_sets: dict[str, list[str]], min_size: int = 5
) -> dict[str, float]:
    """Z-score pathway activity: sum(z_set) / sqrt(|set|).

    If ``z`` is approximately standard normal under the null this is itself a
    z-statistic, so it is comparable between an intervention contrast and a
    patient-vs-reference vector. Sets with fewer than ``min_size`` present genes
    are omitted (not zero-filled).
    """
    idx = {g: i for i, g in enumerate(genes)}
    out: dict[str, float] = {}
    zz = np.nan_to_num(np.asarray(z, dtype=float))
    for name, members in gene_sets.items():
        hit = [idx[g] for g in members if g in idx]
        if len(hit) < min_size:
            continue
        out[name] = float(zz[hit].sum() / np.sqrt(len(hit)))
    return out


def ssgsea(
    genes: list[str], expr: np.ndarray, gene_sets: dict[str, list[str]],
    alpha: float = 0.25, min_size: int = 5,
) -> dict[str, float]:
    """Single-sample GSEA (Barbie et al. 2009) for one expression vector.

    For raw single-sample profiles without a reference cohort. Returns the
    un-normalised enrichment score (sum of the running-sum difference).
    """
    x = np.asarray(expr, dtype=float)
    order = np.argsort(-x, kind="mergesort")
    ranked = [genes[i] for i in order]
    ranks = np.arange(len(x), 0, -1, dtype=float)  # highest expression -> largest rank weight
    n = len(ranked)
    out: dict[str, float] = {}
    for name, members in gene_sets.items():
        mset = set(members)
        in_set = np.array([g in mset for g in ranked])
        k = int(in_set.sum())
        if k < min_size:
            continue
        w = np.where(in_set, ranks**alpha, 0.0)
        p_hit = np.cumsum(w) / w.sum()
        p_miss = np.cumsum(~in_set) / (n - k)
        out[name] = float((p_hit - p_miss).sum())
    return out
