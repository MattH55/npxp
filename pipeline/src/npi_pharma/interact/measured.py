"""Measured interaction signatures from factorial GEO designs.

Every other interaction score in this package is *inferred* from single-agent
signatures, and that approach was measured not to work: signature composition
carries no detectable information about synergy across 50,000 drug-pair
observations (docs/validation_drugcomb.md). A series that profiles control,
agent A, agent B and the A+B combination is different in kind, because the
interaction can be computed directly:

    I(A, B) = (combo - control) - [(A - control) + (B - control)]
            = combo - A - B + control

on a log-expression scale, which is the standard two-factor interaction contrast
and the transcriptional analogue of a Bliss excess. Genes with a large positive
I are induced by the combination beyond what the two agents do separately;
genes with large negative I are suppressed beyond additivity.

This is observed, not assumed. Its limits are the usual ones for a contrast built
from four small groups: it is noisier than any single main effect, since four
group means enter it, and it describes one cell line at one dose and time.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from ..model import Signature
from ..signatures.build import _moderated, standardize

KIND_INTERACTION = "interaction"


def interaction_contrast(control: pd.DataFrame, agent_a: pd.DataFrame, agent_b: pd.DataFrame,
                         combo: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """(effect, moderated statistic) for ``combo - a - b + control`` on shared genes.

    All four frames are genes x samples on a log scale. The standard error combines
    all four groups, so the statistic is shrunk toward zero where any arm is noisy.
    """
    frames = [control, agent_a, agent_b, combo]
    genes = frames[0].index
    for f in frames[1:]:
        genes = genes.intersection(f.index)
    if len(genes) < 200:
        raise ValueError(f"only {len(genes)} genes shared across the four arms")
    arrays = [f.loc[genes].to_numpy(float) for f in frames]
    ns = [a.shape[1] for a in arrays]
    if min(ns) < 2:
        raise ValueError(f"every arm needs >= 2 samples, got {ns}")
    means = [np.nanmean(a, axis=1) for a in arrays]
    variances = [np.nanvar(a, axis=1, ddof=1) for a in arrays]
    effect = means[3] - means[1] - means[2] + means[0]
    se = np.sqrt(sum(v / n for v, n in zip(variances, ns)))
    n_eff = 1.0 / sum(1.0 / n for n in ns)
    stat = _moderated(effect, se * np.sqrt(n_eff), n_eff)
    return pd.Series(effect, index=genes), pd.Series(stat, index=genes)


def build_interaction_signature(
    control: pd.DataFrame, agent_a: pd.DataFrame, agent_b: pd.DataFrame, combo: pd.DataFrame,
    metadata: dict[str, Any],
) -> Signature:
    """A :class:`Signature` of the measured interaction between two agents.

    ``metadata`` should carry ``sig_id``, ``agent_a``, ``agent_b``, ``cell_line`` and
    ``source_accessions``. The stored z is the standardised interaction statistic, so
    it is on the same scale as the drug and NPI signatures and can be compared with
    them; ``meta["effect"]`` keeps the raw log-scale contrast.
    """
    effect, stat = interaction_contrast(control, agent_a, agent_b, combo)
    z = standardize(stat)
    ns = {"n_control": control.shape[1], "n_a": agent_a.shape[1],
          "n_b": agent_b.shape[1], "n_combo": combo.shape[1]}
    flag = "ok" if min(ns.values()) >= 3 else "small_n"
    return Signature(
        sig_id=metadata["sig_id"], kind=KIND_INTERACTION, genes=list(z.index), z=z.to_numpy(),
        provenance=metadata.get("provenance", "geo_factorial"),
        modality="drug_interaction", tissue=metadata.get("cell_line"),
        species=metadata.get("species", "human"),
        sample_size=int(min(ns.values())),
        contrast="combo - a - b + control (log scale), unpaired",
        quality_flag=flag,
        source_accessions=list(metadata.get("source_accessions", [])),
        meta={k: v for k, v in metadata.items() if k not in ("sig_id",)} | ns |
             {"effect": {g: float(v) for g, v in effect.nlargest(50).items()},
              "effect_negative": {g: float(v) for g, v in effect.nsmallest(50).items()}},
    )


def additivity_summary(effect: pd.Series, top: int = 20) -> dict[str, Any]:
    """How far from additive the combination is, and which genes carry it."""
    e = effect.dropna()
    return {
        "n_genes": int(len(e)),
        "sd_of_interaction": float(e.std()),
        "genes_above_additive": {g: round(float(v), 3) for g, v in e.nlargest(top).items()},
        "genes_below_additive": {g: round(float(v), 3) for g, v in e.nsmallest(top).items()},
    }
