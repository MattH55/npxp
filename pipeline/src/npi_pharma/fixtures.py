"""Deterministic SYNTHETIC fixtures for tests and the offline demo.

Nothing here is data. Every signature is tagged provenance="synthetic_fixture"
and every id starts with "FIXTURE_", so fixture output cannot be mistaken for
a real LINCS / GEO result. Signatures are generated from latent pathway
programs over the bundled gene sets so pathway-level logic has something
coherent to find.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .gene_sets import load_gmt
from .signatures.build import build_signature
from .model import DRUG, NPI, PROV_FIXTURE, Signature

SEED = 20260927
N_FILLER = 2400
N_RANDOM_DRUGS = 47

# Program effects (pathway -> signed magnitude). Loosely inspired by known
# biology so the demo reads sensibly; the numbers themselves are invented.
DRUG_PROGRAMS = {
    "FIXTURE_metformin": {"OXIDATIVE_PHOSPHORYLATION_CORE": -2.5, "AMPK_SIGNALING_CORE": 2.0,
                          "MTORC1_SIGNALING_CORE": -1.0, "LIPOGENESIS_CORE": -1.0},
    "FIXTURE_rapamycin": {"MTORC1_SIGNALING_CORE": -3.0, "CELL_CYCLE_E2F_CORE": -1.5,
                          "INSULIN_SIGNALING_CORE": -1.2, "LIPOGENESIS_CORE": -1.0},
    "FIXTURE_pxr_agonist": {"XENOBIOTIC_METABOLISM_CORE": 3.0, "CHOLESTEROL_BIOSYNTHESIS_CORE": 0.5},
}
NPI_PROGRAMS = {
    "FIXTURE_LCD_adipose": ({"LIPOGENESIS_CORE": -2.5, "INFLAMMATORY_RESPONSE_CORE": -1.5,
                             "INSULIN_SIGNALING_CORE": 1.5, "CHOLESTEROL_BIOSYNTHESIS_CORE": -1.0},
                            dict(modality="diet", tissue="adipose", duration="synthetic", intensity="LCD")),
    "FIXTURE_xenobiotic_diet_liver": ({"XENOBIOTIC_METABOLISM_CORE": 2.5, "INFLAMMATORY_RESPONSE_CORE": -0.5},
                                      dict(modality="diet", tissue="liver", duration="synthetic", intensity=None)),
    "FIXTURE_psychosocial_blood": ({"INFLAMMATORY_RESPONSE_CORE": -1.0},
                                   dict(modality="psychosocial", tissue="blood", duration="synthetic",
                                        intensity="synthetic", sample_size=4)),
}
# Built through build_signature from pre/post matrices (same "study", two arms).
STUDY_ARMS = {
    "FIXTURE_CR_muscle": {"MTORC1_SIGNALING_CORE": -2.0, "LIPOGENESIS_CORE": -1.5, "AMPK_SIGNALING_CORE": 1.0},
    "FIXTURE_EX_muscle": {"OXIDATIVE_PHOSPHORYLATION_CORE": 2.5, "FATTY_ACID_OXIDATION_CORE": 2.0,
                          "HYPOXIA_CORE": 1.0, "INFLAMMATORY_RESPONSE_CORE": -0.5},
}
ADIPOSE_DISEASE = {"LIPOGENESIS_CORE": 1.5, "INFLAMMATORY_RESPONSE_CORE": 2.0,
                   "INSULIN_SIGNALING_CORE": -1.5, "OXIDATIVE_PHOSPHORYLATION_CORE": -1.0,
                   "MTORC1_SIGNALING_CORE": 1.0}


class _World:
    def __init__(self, seed: int = SEED):
        self.rng = np.random.default_rng(seed)
        self.sets = load_gmt()
        pathway_genes = sorted({g for m in self.sets.values() for g in m})
        self.genes = pathway_genes + [f"FXG{i:05d}" for i in range(N_FILLER)]
        self.idx = {g: i for i, g in enumerate(self.genes)}
        # Fixed per-gene loadings per program keep pathways coherent across signatures.
        self.loadings = {
            k: np.array([self.idx[g] for g in m]) for k, m in self.sets.items()
        }
        self.weights = {k: self.rng.uniform(0.5, 1.5, len(v)) for k, v in self.loadings.items()}
        self.baseline = self.rng.normal(6.0, 2.0, len(self.genes))

    def program(self, effects: dict[str, float], noise: float = 1.0) -> np.ndarray:
        v = self.rng.normal(0, noise, len(self.genes))
        for k, e in effects.items():
            v[self.loadings[k]] += e * self.weights[k]
        return v


def _sig(world: _World, sid: str, kind: str, z: np.ndarray, **kw) -> Signature:
    return Signature(sig_id=sid, kind=kind, genes=list(world.genes), z=z, provenance=PROV_FIXTURE,
                     quality_flag=kw.pop("quality_flag", "synthetic_fixture"),
                     source_accessions=["SYNTHETIC"], **kw)


def make_drugs(world: _World) -> list[Signature]:
    out = [
        _sig(world, sid, DRUG, world.program(eff), modality="compound", tissue="cell_line_consensus",
             contrast="synthetic", meta={"fixture": True})
        for sid, eff in DRUG_PROGRAMS.items()
    ]
    names = list(world.sets)
    for i in range(N_RANDOM_DRUGS):
        k = world.rng.integers(1, 4)
        eff = {names[j]: float(world.rng.normal(0, 1.8)) for j in world.rng.choice(len(names), k, replace=False)}
        out.append(_sig(world, f"FIXTURE_cmpd_{i:03d}", DRUG, world.program(eff), modality="compound",
                        tissue="cell_line_consensus", contrast="synthetic", meta={"fixture": True, "programs": eff}))
    return out


def make_study_arms(world: _World, n_subjects: int = 10) -> list[Signature]:
    """Mouse CR vs exercise arms of one synthetic study, via the real builder.
    Gene symbols are written mouse-style (Mtor) to exercise normalisation."""
    mouse = [g.capitalize() if not g.startswith("FXG") else g for g in world.genes]
    out = []
    for sid, eff in STUDY_ARMS.items():
        subj_base = world.baseline[:, None] + world.rng.normal(0, 0.8, (len(world.genes), n_subjects))
        pre = subj_base + world.rng.normal(0, 0.3, subj_base.shape)
        shift = np.zeros(len(world.genes))
        for k, e in eff.items():
            shift[world.loadings[k]] += 0.4 * e * world.weights[k]
        post = subj_base + shift[:, None] + world.rng.normal(0, 0.3, subj_base.shape)
        cols = [f"S{i}" for i in range(n_subjects)]
        md = dict(npi_id=sid, modality={"CR": "caloric_restriction", "EX": "exercise"}[sid.split("_")[1]],
                  tissue="skeletal_muscle", species="mouse", duration="synthetic", intensity="synthetic",
                  provenance=PROV_FIXTURE, quality_flag="synthetic_fixture",
                  source_accessions=["SYNTHETIC_CR_EX_STUDY"])
        out.append(build_signature(pd.DataFrame(pre, index=mouse, columns=cols),
                                   pd.DataFrame(post, index=mouse, columns=cols), md))
    return out


def make_npis(world: _World) -> list[Signature]:
    out = []
    for sid, (eff, md) in NPI_PROGRAMS.items():
        md = dict(md)
        n = md.pop("sample_size", 20)
        out.append(_sig(world, sid, NPI, world.program(eff), species="human", sample_size=n,
                        contrast="synthetic", quality_flag="small_n" if n < 6 else "synthetic_fixture", **md))
    return out + make_study_arms(world)


def make_patient_matrices(world: _World, disease: dict[str, float], n_ref: int = 12) -> tuple[pd.DataFrame, pd.DataFrame]:
    ref = world.baseline[:, None] + world.rng.normal(0, 0.5, (len(world.genes), n_ref))
    pat = world.baseline + world.rng.normal(0, 0.5, len(world.genes))
    for k, e in disease.items():
        pat[world.loadings[k]] += 0.5 * e * world.weights[k]
    return (pd.DataFrame({"FIXTURE_patient": pat}, index=world.genes),
            pd.DataFrame(ref, index=world.genes, columns=[f"REF{i:02d}" for i in range(n_ref)]))


def build_all(seed: int = SEED) -> dict[str, object]:
    world = _World(seed)
    npis = make_npis(world)
    drugs = make_drugs(world)
    patient, reference = make_patient_matrices(world, ADIPOSE_DISEASE)
    return {"world": world, "npis": npis, "drugs": drugs, "patient_expr": patient, "reference_expr": reference}


def write_fixture(out_dir: str | Path, seed: int = SEED) -> dict[str, Path]:
    from .store import save_signatures

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fx = build_all(seed)
    paths = {
        "npis": save_signatures(fx["npis"], out / "npis.parquet"),
        "drugs": save_signatures(fx["drugs"], out / "drugs.parquet"),
        "patient": out / "patient_adipose.tsv",
        "reference": out / "reference_adipose.tsv",
    }
    fx["patient_expr"].to_csv(paths["patient"], sep="\t", float_format="%.5f")
    fx["reference_expr"].to_csv(paths["reference"], sep="\t", float_format="%.5f")
    return paths
