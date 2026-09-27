"""Encode an individual transcriptome into a PatientState."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..genes import normalize_symbols
from ..ingest.geo import maybe_log2
from ..model import PatientState


def _clean(expr: pd.DataFrame) -> pd.DataFrame:
    e = expr.copy()
    e.index = normalize_symbols(e.index)
    e = e[~e.index.isin(["", "NAN", "NONE"])]
    e = e.apply(pd.to_numeric, errors="coerce")
    return e.groupby(level=0).mean()


def read_expression(path: str | Path, obs_filter: dict | None = None, sample_key: str | None = None) -> tuple[pd.DataFrame, list[str]]:
    """Read bulk TSV/CSV (genes x samples) or AnnData .h5ad (cells/samples x genes).

    For .h5ad with ``sample_key`` the cells are pseudobulked (mean) per sample;
    ``obs_filter`` (e.g. {"compartment": "malignant"}) selects cells first, the
    hook for scTherapy-style malignant/normal splits.
    """
    path = Path(path)
    flags: list[str] = []
    if path.suffix == ".h5ad":
        import anndata

        ad = anndata.read_h5ad(path)
        for k, v in (obs_filter or {}).items():
            ad = ad[ad.obs[k].astype(str) == str(v)]
        X = ad.X.toarray() if hasattr(ad.X, "toarray") else np.asarray(ad.X)
        df = pd.DataFrame(X, index=ad.obs_names, columns=ad.var_names)
        if sample_key:
            df = df.groupby(ad.obs[sample_key].astype(str).values).mean()
            flags.append("pseudobulk_from_single_cell")
        return df.T, flags
    sep = "," if path.suffix == ".csv" else "\t"
    return pd.read_csv(path, sep=sep, index_col=0), flags


def encode_patient(
    expr: pd.DataFrame,
    sample_id: str | None,
    tissue: str,
    reference: pd.DataFrame | None = None,
    allow_no_reference: bool = False,
    flags: list[str] | None = None,
) -> PatientState:
    """Build s_P for one sample.

    Priority:
      1. ``reference`` = healthy samples of the same tissue/platform: s_P is the
         patient's per-gene z vs that set ("patient - healthy").
      2. >= 3 samples in ``expr``: z vs the cohort (s_P = deviation from cohort,
         not from health; flagged ``cohort_relative``).
      3. single sample, ``allow_no_reference``: within-sample robust z. Mostly
         reflects baseline expression level; flagged and low-confidence.
    """
    flags = list(flags or [])
    e = maybe_log2(_clean(expr))
    if sample_id is None:
        if e.shape[1] != 1:
            raise ValueError(f"expression has {e.shape[1]} samples; pass sample_id")
        sample_id = str(e.columns[0])
    if sample_id not in e.columns:
        raise KeyError(f"sample {sample_id!r} not in expression columns")
    x = e[sample_id]

    if reference is not None:
        r = maybe_log2(_clean(reference))
        genes = x.index.intersection(r.index)
        mu, sd = r.loc[genes].mean(axis=1), r.loc[genes].std(axis=1, ddof=1)
        desc = f"healthy reference (n={r.shape[1]})"
        if r.shape[1] < 5:
            flags.append("small_reference")
    elif e.shape[1] >= 3:
        genes = x.index
        cohort = e.drop(columns=[sample_id]) if e.shape[1] > 3 else e
        mu, sd = cohort.mean(axis=1), cohort.std(axis=1, ddof=1)
        desc = f"cohort (n={cohort.shape[1]}); deviation from cohort, not from health"
        flags.append("cohort_relative")
    elif allow_no_reference:
        genes = x.index
        med = x.median()
        mu = pd.Series(med, index=genes)
        sd = pd.Series(1.4826 * (x - med).abs().median(), index=genes)
        desc = "within-sample robust z (no reference)"
        flags.append("no_reference")
    else:
        raise ValueError("need a healthy reference or >= 3 cohort samples (or allow_no_reference)")

    sd = sd.loc[genes]
    floor = float(np.nanmedian(sd[sd > 0])) if (sd > 0).any() else 1.0
    z = (x.loc[genes] - mu.loc[genes]) / np.maximum(sd, 0.1 * floor)
    z = z.replace([np.inf, -np.inf], np.nan).dropna()
    return PatientState(
        sample_id=str(sample_id),
        tissue=tissue,
        genes=list(z.index),
        expr_z=z.to_numpy(),
        disease_vector=z.to_numpy(),
        reference=desc,
        flags=flags,
    )
