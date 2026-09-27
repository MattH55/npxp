"""Shared gene space: symbol normalisation and alignment across matrices."""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np

MIN_GENE_OVERLAP = 200


class GeneOverlapError(ValueError):
    """Raised when aligned gene spaces are too small to score."""


def normalize_symbol(symbol: str) -> str:
    """Upper-case, strip whitespace and version suffixes.

    Upper-casing also maps mouse symbols (``Cyp3a11`` -> ``CYP3A11``) onto the
    human namespace for the ~80% of 1:1 orthologs that share a symbol; callers
    must flag cross-species comparisons themselves.
    """
    s = str(symbol).strip().upper()
    if "." in s and s.split(".")[-1].isdigit():
        s = s.rsplit(".", 1)[0]
    return s


def normalize_symbols(symbols: Iterable[str]) -> list[str]:
    return [normalize_symbol(s) for s in symbols]


def dedupe_vector(genes: Sequence[str], values: np.ndarray) -> tuple[list[str], np.ndarray]:
    """Normalise symbols and collapse duplicates by keeping the max-|value| entry."""
    norm = normalize_symbols(genes)
    best: dict[str, int] = {}
    for i, g in enumerate(norm):
        if not g or g in {"NAN", "NONE", "---", "-"}:
            continue
        v = values[i]
        if np.isnan(v):
            continue
        j = best.get(g)
        if j is None or abs(v) > abs(values[j]):
            best[g] = i
    out_genes = sorted(best)
    return out_genes, np.array([values[best[g]] for g in out_genes], dtype=float)


def intersect(*gene_lists: Sequence[str]) -> list[str]:
    if not gene_lists:
        return []
    common = set(gene_lists[0])
    for gl in gene_lists[1:]:
        common &= set(gl)
    return sorted(common)


def reindex(genes: Sequence[str], values: np.ndarray, target: Sequence[str]) -> np.ndarray:
    idx = {g: i for i, g in enumerate(genes)}
    return np.array([values[idx[g]] for g in target], dtype=float)


def align(
    named: dict[str, tuple[Sequence[str], np.ndarray]],
    min_overlap: int = MIN_GENE_OVERLAP,
) -> tuple[list[str], dict[str, np.ndarray]]:
    """Align several (genes, vector) pairs onto their intersection.

    Raises :class:`GeneOverlapError` instead of silently returning a tiny or
    empty space.
    """
    shared = intersect(*(g for g, _ in named.values()))
    if len(shared) < min_overlap:
        sizes = ", ".join(f"{k}={len(g)}" for k, (g, _) in named.items())
        raise GeneOverlapError(
            f"shared gene space has {len(shared)} genes (< {min_overlap}); inputs: {sizes}. "
            "Check symbol namespaces (HGNC symbols expected) and platform annotation."
        )
    return shared, {k: reindex(g, v, shared) for k, (g, v) in named.items()}
