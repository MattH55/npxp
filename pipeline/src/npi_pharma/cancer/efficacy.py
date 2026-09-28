"""Combination efficacy from single-agent responses.

Validated on DrugComb's 739,964 measured drug-pair experiments: the Bliss
independence expectation built from each agent's own relative inhibition ranks
the combination's measured efficacy (CSS) at Spearman 0.696 overall and 0.752
within cell line, with no fitting (see docs/efficacy_from_monotherapy.md). A
gradient-boosted model on the same single-agent features plus the pair's mean
reaches 0.809 under leave-cell-lines-out, so the parameter-free form gives up
some accuracy for the ability to transfer to agents that have no training data --
which is the case for every NPI.

What this does NOT do: predict synergy. Measured combination efficacy and Bliss
synergy correlate only 0.060 in DrugComb, and synergy stays near-unpredictable
for unseen agents. These functions estimate how much a combination kills, not
whether it exceeds additivity.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def bliss_expected(a: float | np.ndarray, b: float | np.ndarray) -> float | np.ndarray:
    """Independent (Bliss) expectation of combined inhibition from two single-agent
    inhibitions on a 0-1 scale: ``1 - (1-a)(1-b)``.

    Values are clipped into [0, 1] first: a negative inhibition (an agent that
    increased growth) carries no information about combined kill in this model,
    and a value above 1 is not meaningful.
    """
    a = np.clip(np.asarray(a, dtype=float), 0.0, 1.0)
    b = np.clip(np.asarray(b, dtype=float), 0.0, 1.0)
    out = 1.0 - (1.0 - a) * (1.0 - b)
    return float(out) if np.ndim(out) == 0 else out


def excess_over_bliss(combined: float | np.ndarray, a, b):
    """Measured combined inhibition minus the Bliss expectation (0-1 scale).

    This is the synergy quantity. It is defined here for completeness and for
    scoring measured data; it must not be read off a *predicted* ``combined``.
    """
    return np.asarray(combined, dtype=float) - bliss_expected(a, b)


def auc_to_inhibition(auc: float | np.ndarray) -> float | np.ndarray:
    """Dose-response AUC (1 = no effect over the tested range) -> inhibition 0-1.

    Applies to PRISM and GDSC fitted AUC. This is a monotone rescaling, not a
    calibration: it preserves the ranking of drugs within a cell line, which is
    what the efficacy ranking uses.
    """
    return np.clip(1.0 - np.asarray(auc, dtype=float), 0.0, 1.0)


def rank_combinations(npi_inhibition: pd.Series, drug_inhibition: pd.DataFrame) -> pd.DataFrame:
    """Rank NPI x drug per cell line by expected combined inhibition.

    ``npi_inhibition`` is indexed by cell line id; ``drug_inhibition`` is cell
    lines x drugs. Only cell lines present in both are returned, so a missing NPI
    monotherapy measurement drops that cell line rather than being imputed.
    """
    lines = npi_inhibition.index.intersection(drug_inhibition.index)
    rows = []
    for cl in lines:
        for drug in drug_inhibition.columns:
            di = drug_inhibition.at[cl, drug]
            if pd.isna(di):
                continue
            ni = float(npi_inhibition[cl])
            rows.append({"cell_line_id": cl, "drug": drug, "npi_inhibition": ni,
                         "drug_inhibition": float(di),
                         "expected_combined_inhibition": bliss_expected(ni, di)})
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values("expected_combined_inhibition", ascending=False).reset_index(drop=True)
