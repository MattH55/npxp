"""Label-free NPI x drug interaction scores for one patient (v1 core).

Notation: s_P patient disease vector (what should be reversed), s_N NPI
signature, s_D drug signature, all aligned to one shared gene space.

Outputs are Tier 3 (mechanism-only) plausibility scores: a ranking device,
never a measured or calibrated synergy.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from .. import EVIDENCE_TIER
from ..config import canonical_tissue
from ..gene_sets import pathway_activity
from ..genes import align
from ..model import PROV_FIXTURE, PatientState, Signature

COMPONENTS = ("complementarity", "orthogonality", "monotherapy_correlation", "pathway_joint", "mono_reversal")


def _unit(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(a @ b / (na * nb)) if na > 0 and nb > 0 else 0.0


def reverse(s_p: np.ndarray, s: np.ndarray) -> float:
    """Reversal of the patient state by signature s: -cos(s_P, s) in [-1, 1]."""
    return -cosine(s_p, s)


def combine(s_n: np.ndarray, s_d: np.ndarray, mode: str = "unit_sum") -> np.ndarray:
    if mode == "unit_sum":
        return _unit(s_n) + _unit(s_d)
    if mode == "sum":
        return s_n + s_d
    raise ValueError(f"unknown combine mode {mode!r}")


def pathway_joint(
    a_p: dict[str, float], a_n: dict[str, float], a_d: dict[str, float],
    scale: float, penalty: float,
) -> tuple[float, list[dict[str, Any]]]:
    """DIPx-like pathway compatibility.

    need_k = -tanh(a_P,k / c): direction the patient needs pathway k moved.
    joint_k = need_k * (u_N + u_D)/2 - lambda * |need_k| * conflict_k, with
    u = tanh(a / c) and conflict_k = min(|u_N|, |u_D|) when NPI and drug move k
    in opposite directions. Score = sum(joint) / sum(|need|), clipped to [-1, 1].
    """
    rows, tot, norm = [], 0.0, 0.0
    for k in sorted(set(a_p) & set(a_n) & set(a_d)):
        need = -np.tanh(a_p[k] / scale)
        un, ud = np.tanh(a_n[k] / scale), np.tanh(a_d[k] / scale)
        conflict = min(abs(un), abs(ud)) if un * ud < 0 else 0.0
        j = need * (un + ud) / 2 - penalty * abs(need) * conflict
        tot += j
        norm += abs(need)
        rows.append({
            "pathway": k, "joint": float(j), "patient_activity": a_p[k],
            "npi_activity": a_n[k], "drug_activity": a_d[k], "conflict": float(conflict),
        })
    rows.sort(key=lambda r: -abs(r["joint"]))
    return (float(np.clip(tot / norm, -1, 1)) if norm > 0 else 0.0), rows


def score_pair(
    patient: PatientState,
    npi: Signature,
    drug: Signature,
    gene_sets: dict[str, list[str]],
    cfg: dict[str, Any],
) -> dict[str, Any]:
    genes, v = align(
        {"patient": (patient.genes, patient.disease_vector), "npi": (npi.genes, npi.z), "drug": (drug.genes, drug.z)},
        min_overlap=cfg.get("min_overlap_genes", 200),
    )
    p, n, d = v["patient"], v["npi"], v["drug"]
    combo = combine(n, d, cfg.get("combine", "unit_sum"))

    k = min(cfg.get("focus_top_k", 500), len(genes))
    focus = np.argsort(-np.abs(p))[:k]
    # Reversal is scored on the patient's focus genes (CMap-style query) by
    # default: on the full space it is dominated by genes the patient never moved.
    rs = focus if cfg.get("reversal_space", "focus") == "focus" else slice(None)
    rev_n, rev_d, rev_c = reverse(p[rs], n[rs]), reverse(p[rs], d[rs]), reverse(p[rs], combo[rs])
    best = max(rev_n, rev_d)
    gain = rev_c - best
    # Normalised gain: share of the remaining reversal headroom (1 - best) the
    # combination captures; losses are scaled by the room to fall (1 + best).
    denom = (1.0 - best) if gain >= 0 else (1.0 + best)
    complementarity = float(np.clip(gain / denom, -1, 1)) if denom > 1e-12 else 0.0
    cos_nd = cosine(n, d)
    orth = 1.0 - abs(cosine(n[focus], d[focus]))

    ms = cfg.get("pathway_min_size", 5)
    a_p = pathway_activity(genes, p, gene_sets, ms)
    a_n = pathway_activity(genes, n, gene_sets, ms)
    a_d = pathway_activity(genes, d, gene_sets, ms)
    pj, pw_rows = pathway_joint(a_p, a_n, a_d, cfg.get("pathway_scale", 3.0), cfg.get("pathway_conflict_penalty", 1.0))

    scores = {
        "complementarity": complementarity,
        # Orthogonality only counts toward the composite when both agents reverse
        # s_P ("orthogonal and both helpful"); unrelated drugs are trivially orthogonal.
        "orthogonality": orth if (rev_n > 0 and rev_d > 0) or not cfg.get("gate_orthogonality", True) else 0.0,
        "monotherapy_correlation": cos_nd,
        "pathway_joint": pj,
        "mono_reversal": (rev_n + rev_d) / 2,
    }
    w = cfg["weights"]
    wsum = sum(abs(w.get(c, 0.0)) for c in COMPONENTS) or 1.0
    composite_raw = sum(w.get(c, 0.0) * scores[c] for c in COMPONENTS) / wsum

    flags, warnings, confidence = [], [], 1.0
    shrink = cfg.get("shrink", {})

    def penalize(key: str, msg: str) -> None:
        nonlocal confidence
        flags.append(key)
        warnings.append(msg)
        confidence *= shrink.get(key, 1.0)

    pt, nt = canonical_tissue(patient.tissue, cfg), canonical_tissue(npi.tissue, cfg)
    if pt != nt:
        penalize("tissue_mismatch", f"patient tissue '{patient.tissue}' vs NPI signature tissue '{npi.tissue}'; "
                 "NPI programs are tissue-specific, score down-weighted")
    if npi.quality_flag in ("small_n", "exploratory", "unverified_metadata"):
        penalize(npi.quality_flag, f"NPI signature quality_flag={npi.quality_flag} (n={npi.sample_size})")
    if npi.species != "human" or drug.species != "human":
        penalize("cross_species", f"cross-species comparison (NPI {npi.species}, drug {drug.species}); "
                 "symbols matched by upper-casing, not ortholog mapping")
    if len(genes) < cfg.get("low_overlap_genes", 1000):
        penalize("low_overlap", f"only {len(genes)} shared genes")
    if PROV_FIXTURE in (npi.provenance, drug.provenance):
        flags.append("SYNTHETIC_FIXTURE")
        warnings.append("one or both signatures are synthetic test fixtures, not data")

    fc = cfg.get("flags", {})
    thr = fc.get("pathway_activity", 2.0)
    if cos_nd > fc.get("redundancy_cos", 0.6) and rev_n > 0 and rev_d > 0:
        flags.append("possible_redundancy")
        warnings.append(f"NPI and drug programs highly correlated (cos={cos_nd:.2f}); may be redundant rather than additive")
    antagonism = []
    for r in pw_rows:
        if (abs(r["patient_activity"]) > thr and abs(r["npi_activity"]) > thr and abs(r["drug_activity"]) > thr
                and np.sign(r["npi_activity"]) != np.sign(r["drug_activity"])):
            antagonism.append({
                "pathway": r["pathway"],
                "patient_needs": "down" if r["patient_activity"] > 0 else "up",
                "npi": "up" if r["npi_activity"] > 0 else "down",
                "drug": "up" if r["drug_activity"] > 0 else "down",
            })
    if antagonism:
        flags.append("pathway_antagonism")
    xs = fc.get("xenobiotic_set", "XENOBIOTIC_METABOLISM_CORE")
    if a_n.get(xs, 0.0) > thr and abs(a_d.get(xs, 0.0)) > thr:
        flags.append("xenobiotic_genes_moved_by_both")
        warnings.append(f"NPI induces {xs} genes and the drug also moves them: a possible metabolic (PK) "
                        "interaction worth checking. Transcriptomic observation only, not a PK prediction.")

    # Gene attributions: per-gene terms of rev(combo) = sum_g -p_hat_g * c_hat_g.
    pr, cr = np.zeros(len(genes)), np.zeros(len(genes))
    pr[rs], cr[rs] = p[rs], combo[rs]
    contrib = -_unit(pr) * _unit(cr)
    order = np.argsort(-contrib)
    na = cfg.get("n_attribution_genes", 20)

    def gene_row(i: int) -> dict[str, Any]:
        return {"gene": genes[i], "contribution": float(contrib[i]), "patient_z": float(p[i]),
                "npi_z": float(n[i]), "drug_z": float(d[i])}

    return {
        "patient": patient.sample_id,
        "npi_id": npi.sig_id,
        "drug_id": drug.sig_id,
        "evidence_tier": EVIDENCE_TIER,
        "n_genes": len(genes),
        **scores,
        "complementarity_gain": float(gain),
        "orthogonality_raw": orth,
        "reverse_npi": rev_n,
        "reverse_drug": rev_d,
        "reverse_combo": rev_c,
        "composite_raw": float(composite_raw),
        "confidence": float(confidence),
        "composite": float(composite_raw * confidence),
        "flags": flags + [f"patient:{f}" for f in patient.flags],
        "warnings": warnings,
        "antagonism": antagonism,
        "pathways": pw_rows,
        "genes_driving_reversal": [gene_row(i) for i in order[:na]],
        "genes_opposing_reversal": [gene_row(i) for i in order[::-1][:na]],
        "weights": {c: w.get(c, 0.0) for c in COMPONENTS},
    }
