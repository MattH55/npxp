"""Core data objects: intervention signatures and patient states."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

NPI = "npi"
DRUG = "drug"

# Provenance values. Anything tagged synthetic_fixture must never be read as data.
PROV_LINCS = "lincs_l1000"
PROV_GEO = "geo"
PROV_LOCAL = "local_matrix"
PROV_FIXTURE = "synthetic_fixture"


@dataclass
class Signature:
    sig_id: str
    kind: str  # "npi" | "drug"
    genes: list[str]
    z: np.ndarray
    provenance: str
    modality: str | None = None  # NPI: diet/exercise/CR/...; drug: "compound"
    tissue: str | None = None  # NPI tissue, or drug tissue proxy (cell-line lineage)
    species: str = "human"
    duration: str | None = None
    intensity: str | None = None
    sample_size: int | None = None
    contrast: str | None = None
    quality_flag: str = "ok"  # ok | small_n | exploratory | unverified_metadata | synthetic_fixture
    source_accessions: list[str] = field(default_factory=list)
    up_genes: list[str] = field(default_factory=list)
    down_genes: list[str] = field(default_factory=list)
    meta: dict[str, Any] = field(default_factory=dict)  # drug: pert_id, targets, cell_line, dose, time_h, tas, smiles

    def __post_init__(self) -> None:
        self.z = np.asarray(self.z, dtype=float)
        if len(self.genes) != len(self.z):
            raise ValueError(f"{self.sig_id}: {len(self.genes)} genes vs {len(self.z)} values")
        if not self.up_genes and not self.down_genes and len(self.z):
            self.up_genes, self.down_genes = top_genes(self.genes, self.z)

    def as_series(self):
        import pandas as pd

        return pd.Series(self.z, index=self.genes, name=self.sig_id)


def top_genes(genes: list[str], z: np.ndarray, n: int = 150) -> tuple[list[str], list[str]]:
    order = np.argsort(z)
    n = min(n, len(z) // 2)
    up = [genes[i] for i in order[::-1][:n] if z[i] > 0]
    down = [genes[i] for i in order[:n] if z[i] < 0]
    return up, down


@dataclass
class PatientState:
    sample_id: str
    tissue: str
    genes: list[str]
    expr_z: np.ndarray  # expression z-scored vs reference/cohort
    disease_vector: np.ndarray  # s_P: the state an intervention should reverse
    reference: str  # description of what s_P is relative to
    flags: list[str] = field(default_factory=list)
    meta: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.expr_z = np.asarray(self.expr_z, dtype=float)
        self.disease_vector = np.asarray(self.disease_vector, dtype=float)
