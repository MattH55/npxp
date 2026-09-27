"""GEO series-matrix parsing and pre/post pairing for NPI studies."""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd


def _open(path: str | Path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def read_series_matrix(path: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse a GEO ``*_series_matrix.txt(.gz)``.

    Returns (expression probes x GSM, sample annotations GSM x fields). Each
    ``!Sample_characteristics_ch1`` line of the form ``key: value`` becomes a
    column named ``key``.
    """
    ann: dict[str, list[str]] = {}
    rows: list[list[str]] = []
    header: list[str] | None = None
    in_table = False
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                parts = [p.strip('"') for p in line.split("\t")]
                if header is None:
                    header = parts
                else:
                    rows.append(parts)
                continue
            if line.startswith("!Sample_"):
                parts = [p.strip('"') for p in line.split("\t")]
                key, vals = parts[0][len("!Sample_"):], parts[1:]
                if key.startswith("characteristics"):
                    by_key: dict[str, list[str]] = {}
                    for v in vals:
                        k, _, val = v.partition(":")
                        by_key.setdefault(k.strip(), []).append(val.strip())
                    for k, lst in by_key.items():
                        if len(lst) == len(vals):
                            ann[k] = lst
                else:
                    ann.setdefault(key, vals)
    if header is None:
        raise ValueError(f"{path}: no series_matrix table found")
    expr = pd.DataFrame(rows, columns=header).set_index(header[0])
    expr = expr.apply(pd.to_numeric, errors="coerce")
    samples = pd.DataFrame({k: v for k, v in ann.items() if len(v) == expr.shape[1]})
    samples.index = expr.columns
    return expr, samples


def collapse_probes(expr: pd.DataFrame, probe_map: pd.Series) -> pd.DataFrame:
    """Probes -> gene symbols; multiple probes per gene keep the highest-mean probe."""
    m = probe_map.dropna()
    m = m[m.astype(str).str.len() > 0]
    e = expr.loc[expr.index.intersection(m.index)].copy()
    e["__gene"] = m.loc[e.index].astype(str).str.split(" /// ").str[0]
    e["__mean"] = e.drop(columns="__gene").mean(axis=1)
    e = e.sort_values("__mean", ascending=False).drop_duplicates("__gene")
    return e.set_index("__gene").drop(columns="__mean")


def maybe_log2(expr: pd.DataFrame) -> pd.DataFrame:
    """log2(x+1) if the matrix looks linear-scale (GEO2R's heuristic)."""
    q = np.nanquantile(expr.to_numpy(float), [0.0, 0.25, 0.5, 0.99, 1.0])
    if q[4] > 100 or (q[3] - q[0] > 50 and q[1] > 0):
        return np.log2(expr.clip(lower=0) + 1)
    return expr


def split_pre_post(
    expr: pd.DataFrame, samples: pd.DataFrame, pairing: dict
) -> tuple[pd.DataFrame, pd.DataFrame, bool]:
    """Split into pre/post matrices using a catalog ``pairing`` spec:

    ``{time_field, pre, post, subject_field?, filter?: {field: value}}``.
    With ``subject_field`` the columns are matched per subject (paired);
    otherwise returns an unpaired split.
    """
    s = samples.copy()
    for field, value in (pairing.get("filter") or {}).items():
        s = s[s[field].astype(str) == str(value)]
    tf = pairing["time_field"]
    pre_s = s[s[tf].astype(str) == str(pairing["pre"])]
    post_s = s[s[tf].astype(str) == str(pairing["post"])]
    subj = pairing.get("subject_field")
    if subj:
        pre_by = pre_s.reset_index().drop_duplicates(subj).set_index(subj)["index"]
        post_by = post_s.reset_index().drop_duplicates(subj).set_index(subj)["index"]
        both = sorted(set(pre_by.index) & set(post_by.index))
        pre = expr[pre_by.loc[both].tolist()]
        post = expr[post_by.loc[both].tolist()]
        pre.columns = post.columns = both
        return pre, post, True
    return expr[pre_s.index.tolist()], expr[post_s.index.tolist()], False
