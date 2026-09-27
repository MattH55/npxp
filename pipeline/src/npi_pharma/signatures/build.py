"""Build intervention signatures from pre/post expression."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from ..genes import dedupe_vector
from ..model import NPI, Signature

TOP_N = 150


def _moderated(stat_num: np.ndarray, sd: np.ndarray, n_eff: float) -> np.ndarray:
    """t-like statistic with a variance floor (s0 = median gene SD).

    A light-weight stand-in for limma's empirical-Bayes shrinkage: stops
    near-zero-variance genes from dominating the top of the signature.
    """
    s0 = np.nanmedian(sd[sd > 0]) if np.any(sd > 0) else 1.0
    return stat_num / ((sd + s0) / np.sqrt(n_eff))


def paired_de(pre: pd.DataFrame, post: pd.DataFrame) -> pd.Series:
    """Paired DE. ``pre``/``post`` are genes x subjects with matching columns."""
    if list(pre.columns) != list(post.columns):
        common = [c for c in pre.columns if c in set(post.columns)]
        pre, post = pre[common], post[common]
    common_genes = pre.index.intersection(post.index)
    d = (post.loc[common_genes] - pre.loc[common_genes]).to_numpy(float)
    n = d.shape[1]
    if n < 2:
        raise ValueError("paired DE needs >= 2 paired subjects")
    mean = np.nanmean(d, axis=1)
    sd = np.nanstd(d, axis=1, ddof=1)
    return pd.Series(_moderated(mean, sd, n), index=common_genes)


def cohort_de(pre: pd.DataFrame, post: pd.DataFrame) -> pd.Series:
    """Unpaired Welch-style contrast, used when subjects cannot be matched."""
    common_genes = pre.index.intersection(post.index)
    a = pre.loc[common_genes].to_numpy(float)
    b = post.loc[common_genes].to_numpy(float)
    na, nb = a.shape[1], b.shape[1]
    if na < 2 or nb < 2:
        raise ValueError("cohort DE needs >= 2 samples per group")
    diff = np.nanmean(b, axis=1) - np.nanmean(a, axis=1)
    se = np.sqrt(np.nanvar(a, axis=1, ddof=1) / na + np.nanvar(b, axis=1, ddof=1) / nb)
    n_eff = 1.0 / (1.0 / na + 1.0 / nb)
    return pd.Series(_moderated(diff, se * np.sqrt(n_eff), n_eff), index=common_genes)


def standardize(stat: pd.Series) -> pd.Series:
    """Robust z: centre on median, scale by 1.4826*MAD, so signatures from
    different studies/platforms live on a comparable scale."""
    x = stat.astype(float)
    med = np.nanmedian(x)
    mad = 1.4826 * np.nanmedian(np.abs(x - med))
    return (x - med) / (mad if mad > 0 else (np.nanstd(x) or 1.0))


def build_signature(
    expr_pre: pd.DataFrame,
    expr_post: pd.DataFrame,
    metadata: dict[str, Any],
    paired: bool = True,
    min_n: int = 6,
) -> Signature:
    """Build an NPI :class:`Signature` from genes x samples pre/post matrices.

    ``metadata`` must include ``npi_id`` and should carry the catalog fields
    (modality, tissue, species, duration, intensity, contrast, source_accessions,
    provenance). ``sample_size`` defaults to the number of paired subjects.
    """
    if paired:
        stat = paired_de(expr_pre, expr_post)
        n = int(min(expr_pre.shape[1], expr_post.shape[1]))
        contrast = metadata.get("contrast") or "paired post - pre"
    else:
        stat = cohort_de(expr_pre, expr_post)
        n = int(min(expr_pre.shape[1], expr_post.shape[1]))
        contrast = metadata.get("contrast") or "cohort post - pre (unpaired)"
    z = standardize(stat)
    genes, vals = dedupe_vector(list(z.index), z.to_numpy())

    sample_size = metadata.get("sample_size") or n
    flag = metadata.get("quality_flag") or "ok"
    if flag == "ok" and sample_size < min_n:
        flag = "small_n"
    return Signature(
        sig_id=metadata["npi_id"],
        kind=NPI,
        genes=genes,
        z=vals,
        provenance=metadata.get("provenance", "local_matrix"),
        modality=metadata.get("modality"),
        tissue=metadata.get("tissue"),
        species=metadata.get("species", "human"),
        duration=metadata.get("duration"),
        intensity=metadata.get("intensity"),
        sample_size=int(sample_size),
        contrast=contrast,
        quality_flag=flag,
        source_accessions=list(metadata.get("source_accessions", [])),
        meta={k: v for k, v in metadata.items() if k in ("paper", "notes", "arm", "pairing")},
    )


def consensus(sigs: list[Signature], new_id: str) -> Signature:
    """Median consensus of signatures sharing (modality, tissue). Refuses to mix tissues."""
    keys = {(s.modality, s.tissue, s.kind) for s in sigs}
    if len(keys) != 1:
        raise ValueError(f"refusing consensus across modality/tissue/kind: {sorted(map(str, keys))}")
    frame = pd.concat([s.as_series() for s in sigs], axis=1, join="inner")
    med = frame.median(axis=1)
    first = sigs[0]
    return Signature(
        sig_id=new_id, kind=first.kind, genes=list(med.index), z=med.to_numpy(),
        provenance=first.provenance, modality=first.modality, tissue=first.tissue,
        species=first.species, duration=first.duration, intensity=first.intensity,
        sample_size=sum(s.sample_size or 0 for s in sigs) or None,
        contrast=f"median consensus of {len(sigs)}",
        quality_flag=first.quality_flag if all(s.quality_flag == first.quality_flag for s in sigs) else "mixed",
        source_accessions=sorted({a for s in sigs for a in s.source_accessions}),
        meta={"members": [s.sig_id for s in sigs]},
    )
