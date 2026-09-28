"""Drug resistance programs from baseline expression and measured drug response.

A drug's resistance program is, for each gene, the correlation across cancer
cell lines between the gene's baseline expression and the drug's AUC (higher
AUC = more viable = more resistant). Expression and AUC are both centred within
lineage first, so the program is not simply "which tissue resists this drug".

Sources (all public):
  expression  DepMap 24Q4 OmicsExpressionProteinCodingGenesTPMLogp1.csv (figshare 27993248)
  response    PRISM Repurposing 19Q4 secondary-screen-dose-response-curve-parameters.csv
              (figshare 9393293); GDSC1/GDSC2 fitted dose response release 8.2 (Sanger FTP)
  lineage     DepMap Model.csv (OncotreeLineage; SangerModelID for GDSC)
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from ..genes import normalize_symbols


def load_expression(path: str | Path, genes: list[str] | None = None) -> pd.DataFrame:
    """DepMap expression -> cell lines (ModelID) x gene symbols, optionally restricted to ``genes``."""
    header = pd.read_csv(path, nrows=0).columns
    sym = {c: re.sub(r"\s*\(\d+\)$", "", c) for c in header[1:]}
    keep = [c for c in header[1:] if genes is None or sym[c] in set(genes)]
    df = pd.read_csv(path, usecols=[header[0], *keep], index_col=0)
    df.columns = normalize_symbols([sym[c] for c in df.columns])
    df.index.name = "ModelID"
    return df.loc[:, ~df.columns.duplicated()]


def load_prism_auc(path: str | Path, screen: str = "HTS002") -> pd.DataFrame:
    """PRISM 19Q4 secondary AUC -> cell lines x drug names (lower-case), one screen only."""
    p = pd.read_csv(path, usecols=["depmap_id", "name", "auc", "screen_id"])
    p = p[(p["screen_id"] == screen) & p["auc"].notna() & p["name"].notna()]
    p["name"] = p["name"].str.lower()
    return p.pivot_table(index="depmap_id", columns="name", values="auc", aggfunc="median")


def load_gdsc_auc(path: str | Path, models: pd.DataFrame) -> pd.DataFrame:
    """GDSC fitted dose response -> DepMap ModelID x drug names (lower-case); AUC, median over drug ids."""
    g = pd.read_csv(path, usecols=["DRUG_NAME", "SANGER_MODEL_ID", "AUC"])
    m = models.dropna(subset=["SangerModelID"]).drop_duplicates("SangerModelID")
    to_ach = m.set_index("SangerModelID")["ModelID"]
    g["ModelID"] = g["SANGER_MODEL_ID"].map(to_ach)
    g = g.dropna(subset=["ModelID"])
    g["DRUG_NAME"] = g["DRUG_NAME"].str.lower()
    return g.pivot_table(index="ModelID", columns="DRUG_NAME", values="AUC", aggfunc="median")


def center_within(df: pd.DataFrame, groups: pd.Series, min_group: int = 5) -> pd.DataFrame:
    """Subtract group means; rows in groups smaller than ``min_group`` are dropped."""
    g = groups.reindex(df.index)
    sizes = g.map(g.value_counts())
    keep = g.notna() & (sizes >= min_group)
    df, g = df[keep], g[keep]
    return df - df.groupby(g).transform("mean")


def resistance_programs(
    expr: pd.DataFrame, auc: pd.DataFrame, lineage: pd.Series, min_lines: int = 100
) -> tuple[pd.DataFrame, pd.Series]:
    """Genes x drugs matrix of lineage-adjusted Pearson correlations of expression with AUC.

    Each drug uses the cell lines where it was measured; drugs with fewer than
    ``min_lines`` lineage-adjustable lines are skipped. Returns (programs, n_lines).
    """
    lines = expr.index.intersection(auc.index)
    x_all, a_all = expr.loc[lines], auc.loc[lines]
    progs, ns = {}, {}
    for drug in a_all.columns:
        a = a_all[drug].dropna()
        x = x_all.loc[a.index]
        xc = center_within(x, lineage)
        if len(xc) < min_lines:
            continue
        ac = center_within(a.to_frame(), lineage).iloc[:, 0].loc[xc.index]
        xs = (xc - xc.mean()) / xc.std(ddof=0).replace(0, np.nan)
        as_ = (ac - ac.mean()) / (ac.std(ddof=0) or np.nan)
        progs[drug] = (xs.mul(as_, axis=0)).mean()
        ns[drug] = len(xc)
    return pd.DataFrame(progs), pd.Series(ns, name="n_lines")


def split_half_reliability(
    expr: pd.DataFrame, auc: pd.Series, lineage: pd.Series, n_splits: int = 10, seed: int = 0
) -> float:
    """Median correlation between programs built on random halves of the cell lines."""
    rng = np.random.default_rng(seed)
    a = auc.dropna()
    idx = a.index.intersection(expr.index).to_numpy()
    rs = []
    for _ in range(n_splits):
        perm = rng.permutation(idx)
        h1, h2 = perm[: len(perm) // 2], perm[len(perm) // 2:]
        p1, _ = resistance_programs(expr.loc[h1], a.loc[h1].to_frame("d"), lineage, min_lines=20)
        p2, _ = resistance_programs(expr.loc[h2], a.loc[h2].to_frame("d"), lineage, min_lines=20)
        if "d" in p1 and "d" in p2:
            rs.append(p1["d"].corr(p2["d"]))
    return float(np.median(rs)) if rs else float("nan")


def drug_specific(programs: pd.DataFrame, background: pd.DataFrame | None = None) -> pd.DataFrame:
    """Remove the general-sensitivity axis: subtract each program's projection on the
    mean program across ``background`` drugs (default: all columns of ``programs``)."""
    bg = (background if background is not None else programs).mean(axis=1)
    bg = bg.loc[programs.index]
    u = bg / np.linalg.norm(bg)
    return programs - np.outer(u, u.to_numpy() @ programs.to_numpy())
