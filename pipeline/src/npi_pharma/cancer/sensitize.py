"""Score NPI x drug sensitisation as reversal of the drug's resistance program.

score(N, d) = -corr_genes(s_N, r_d): positive when the NPI signature s_N moves
expression toward the state of cell lines that are sensitive to drug d.
``percentile`` places each drug among all scored drugs for that NPI.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def score_matrix(npi_sigs: pd.DataFrame, programs: pd.DataFrame, min_genes: int = 200) -> pd.DataFrame:
    """NPIs (columns of ``npi_sigs``, genes x NPIs) x drugs (columns of ``programs``)."""
    genes = npi_sigs.index.intersection(programs.index)
    if len(genes) < min_genes:
        raise ValueError(f"only {len(genes)} shared genes (< {min_genes})")
    s = npi_sigs.loc[genes].apply(lambda c: (c - c.mean()) / c.std(ddof=0))
    r = programs.loc[genes].apply(lambda c: (c - c.mean()) / c.std(ddof=0))
    return -(s.T.fillna(0) @ r.fillna(0)) / len(genes)


def percentile(scores: pd.DataFrame) -> pd.DataFrame:
    """Per NPI (row), each drug's percentile among all drugs (1 = most sensitising)."""
    return scores.rank(axis=1, pct=True)
