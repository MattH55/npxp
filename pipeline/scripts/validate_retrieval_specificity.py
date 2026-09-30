"""Is an NPI-drug match specific, or an artefact of one dominant axis?

A retrieval result is normally read down the drug axis: of all the drug signatures,
which one best matches this NPI? That question alone is not enough, and this project
was misled by it. Cystine deprivation in MCF7 ranked erastin first of 1,803 drug
signatures, which is mechanistically seductive -- erastin inhibits SLC7A11 and cystine
deprivation starves that same transporter -- and it does not survive two checks.

**The reciprocal check.** Ask the same question down the NPI axis: of all the NPIs,
which does this drug best match? For a mechanism-specific pairing the answer has to be
symmetric. It is not. The LINCS Phase I erastin consensus (51 cell lines) puts cystine
deprivation **6th of 14** NPIs, and even the GEO erastin consensus puts it **3rd**,
both behind serum-free LoVo and glucose deprivation.

**The shared-axis check.** Correlate the drugs' NPI-profiles with each other. If every
drug ranks the NPIs in nearly the same order, that order belongs to the NPIs and
carries no drug-specific information. Measured: median Spearman +0.75 between drugs in
the Phase I panel, with one axis explaining 80% of the variance. Azacitidine's profile
is nearly indistinguishable from erastin's.

What that axis is: the NPIs differ enormously in how much transcriptional response
they produce. The spread of an NPI's similarity across the 1,763-drug LINCS panel runs
from 0.144 (serum-free, 96 h) down to 0.010 for the three amino-acid deprivations --
fourteenfold. A drug's ranking of NPIs mostly recovers that magnitude ordering.

So a top-ranked pairing needs all three: a high rank down the drug axis, a high rank
down the NPI axis, and a similarity not explained by the dominant axis. This script
reports all three.

    python scripts/validate_retrieval_specificity.py
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def profile_matrix(t: pd.DataFrame, panel: str) -> pd.DataFrame:
    """NPIs x drugs similarity for one panel."""
    return t[t["source"] == panel].pivot_table(index="npi", columns="drug",
                                               values="similarity")


def shared_axis(M: pd.DataFrame) -> dict:
    """How much of the panel is one axis shared by every drug?"""
    C = M.corr(method="spearman")
    iu = np.triu_indices_from(C, 1)
    v = C.to_numpy()[iu]
    v = v[np.isfinite(v)]
    X = M.fillna(0.0).to_numpy()
    X = X - X.mean(axis=0)
    s = np.linalg.svd(X, compute_uv=False) if min(X.shape) > 1 else np.array([1.0])
    return {"n_drugs": int(M.shape[1]), "n_npis": int(M.shape[0]),
            "median_between_drug_profile_spearman": float(np.median(v)) if len(v) else None,
            "first_axis_variance_fraction": float(s[0] ** 2 / (s ** 2).sum())}


def two_way_residual(M: pd.DataFrame) -> pd.DataFrame:
    """Similarity double-standardised, by NPI and then by drug.

    Subtracting row and column *means* is not enough here, and a first attempt at this
    gate that did only that passed 229 pairings -- almost all of them on the two NPIs
    with the largest response, and almost all with "this drug's best NPI" rank 1. The
    dominant axis is NPI response *magnitude*, which is a difference in row variance,
    not row mean: the spread of an NPI's similarity across the LINCS panel runs from
    0.145 down to 0.010. Removing a mean leaves that untouched.

    So each NPI's similarities are z-scored across drugs (which removes its mean AND
    its scale), and those z-scores are then z-scored within each drug across NPIs
    (removing a drug's overall tendency to score high). What survives is specific to
    the pairing: large only if the drug is unusual for that NPI *and* the NPI is
    unusual for that drug, on a scale both can be compared on.
    """
    A = M.copy()
    rsd = A.std(axis=1).replace(0, np.nan)
    Z = A.sub(A.mean(axis=1), axis=0).div(rsd, axis=0)
    csd = Z.std(axis=0).replace(0, np.nan)
    return Z.sub(Z.mean(axis=0), axis=1).div(csd, axis=1)


CALIBRATION_NOTE = (
    "CALIBRATION IS AN OPEN ITEM. residual_z is descriptive, not a significance test. "
    "Two reasons. It saturates: the cell inflates the standard deviation of its own "
    "column, so with n NPIs it cannot exceed about sqrt(n - 1) (3.6 here) however large "
    "the real effect -- planting 6, 12, 20 and 40 sd into a test matrix scores 1.98, "
    "2.22, 2.27 and 2.29. And the obvious null does not work: shuffling drugs within "
    "each NPI preserves that NPI's values, so an extreme cell is still extreme in every "
    "shuffle and inflates the null it is tested against -- a planted 200 sd effect "
    "scores p = 0.37 against it. A usable test needs a statistic that cannot mask "
    "itself (deleted or robust standardisation) with a null that does not carry the "
    "effect. Until then, read the two axis ranks, which need no calibration."
)


def sweep(t: pd.DataFrame, min_drugs: int = 8, top: int = 12) -> pd.DataFrame:
    """Every pairing scored on both axes and on the two-way residual."""
    rows = []
    for panel, g in t.groupby("source"):
        M = g.pivot_table(index="npi", columns="drug", values="similarity")
        if M.shape[1] < min_drugs:
            continue
        R = two_way_residual(M.fillna(M.stack().mean()))
        # Rank down each axis on the row-standardised values, so "which NPI does this
        # drug prefer" is not just "which NPI responds most".
        rsd = M.std(axis=1).replace(0, np.nan)
        Zr = M.sub(M.mean(axis=1), axis=0).div(rsd, axis=0)
        drug_rank = M.rank(axis=1, ascending=False)     # within an NPI, across drugs
        npi_rank = Zr.rank(axis=0, ascending=False)     # within a drug, across NPIs
        for npi in M.index:
            for drug in M.columns:
                if pd.isna(M.loc[npi, drug]):
                    continue
                rows.append({
                    "panel": panel, "npi": npi, "drug": drug,
                    "similarity": float(M.loc[npi, drug]),
                    "rank_drug_axis": int(drug_rank.loc[npi, drug]),
                    "n_drugs": int(M.shape[1]),
                    "rank_npi_axis": int(npi_rank.loc[npi, drug]),
                    "n_npis": int(M.shape[0]),
                    "residual_z": float(R.loc[npi, drug]),
                })
    r = pd.DataFrame(rows)
    if not len(r):
        return r
    # A pairing is specific only if it is top-decile down BOTH axes and its residual
    # survives the main-effect ablation.
    r["top_decile_both_axes"] = (
        (r["rank_drug_axis"] <= np.ceil(r["n_drugs"] * 0.1).clip(lower=1)) &
        (r["rank_npi_axis"] <= np.ceil(r["n_npis"] * 0.1).clip(lower=1)))
    # NOT a significance test -- see CALIBRATION_NOTE. A flag for inspection only.
    r["flag_for_inspection"] = r["top_decile_both_axes"] & (r["residual_z"] >= 3.0)
    return r.sort_values("residual_z", ascending=False)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ranking", default="out/npi_drug_all/npi_drug_all.tsv")
    ap.add_argument("--pairs", nargs="*",
                    default=["gse62673_mcf7_cystine_deprivation24h:erastin"],
                    help="NPI:drug pairings to check, as claimed results")
    ap.add_argument("--out", default="out/retrieval_specificity")
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t = pd.read_csv(a.ranking, sep="\t")

    axes = {p: shared_axis(profile_matrix(t, p)) for p in sorted(t["source"].unique())}
    magnitude = (t[t["reliability"] == "high"].groupby("npi")["similarity"].std()
                 .sort_values(ascending=False))

    rows = []
    for spec in a.pairs:
        npi, drug = spec.split(":", 1)
        for panel, g in t[t["drug"] == drug].groupby("source"):
            # rank down the DRUG axis: this NPI against every drug in the panel
            d_axis = t[(t["npi"] == npi) & (t["source"] == panel)]
            d_axis = d_axis.sort_values("similarity", ascending=False).reset_index(drop=True)
            hit = d_axis.index[d_axis["drug"] == drug]
            # rank down the NPI axis: this drug against every NPI
            n_axis = g.sort_values("similarity", ascending=False).reset_index(drop=True)
            nhit = n_axis.index[n_axis["npi"] == npi]
            rows.append({
                "npi": npi, "drug": drug, "panel": panel,
                "similarity": float(g.loc[g["npi"] == npi, "similarity"].iloc[0])
                if (g["npi"] == npi).any() else None,
                "rank_down_drug_axis": int(hit[0]) + 1 if len(hit) else None,
                "n_drugs_in_panel": len(d_axis),
                "rank_down_npi_axis": int(nhit[0]) + 1 if len(nhit) else None,
                "n_npis": len(n_axis),
                "npi_response_magnitude": float(magnitude.get(npi, float("nan"))),
            })
    r = pd.DataFrame(rows)
    r.to_csv(out / "specificity.tsv", sep="\t", index=False)
    print("Claimed pairings, checked down BOTH axes:\n")
    print(r.to_string(index=False, na_rep="-"))
    print("\nOne shared axis per panel (if drugs agree on the NPI ordering, that "
          "ordering is not about the drugs):")
    print(pd.DataFrame(axes).T.to_string(float_format=lambda v: f"{v:.3f}"))
    print("\nNPI response magnitude (sd of similarity across the LINCS panel):")
    print(magnitude.to_string(float_format=lambda v: f"{v:.4f}"))

    # the retrospective sweep over every pairing
    sw = sweep(t)
    if len(sw):
        sw.to_csv(out / "sweep.tsv", sep="\t", index=False)
        passed = sw[sw["flag_for_inspection"]]
        print(f"\n=== retrospective gate over {len(sw)} pairings in "
              f"{sw['panel'].nunique()} panels ===")
        print("gate: top decile down BOTH axes AND residual z >= 3 after removing the "
              "NPI and drug main effects")
        # With this many comparisons a z threshold is met by chance many times over,
        # so the count only means something against that expectation.
        print(f"\n{len(passed)} pairing(s) flagged for inspection (NOT significant):")
        cols = ["panel", "npi", "drug", "similarity", "rank_drug_axis", "n_drugs",
                "rank_npi_axis", "n_npis", "residual_z"]
        print(passed[cols].to_string(index=False) if len(passed) else "  (none)")
        print("\nhighest residual_z regardless of the axis ranks:")
        print(sw.head(10)[cols + ["top_decile_both_axes"]].to_string(index=False))

    verdict = {
        "pairs": rows, "shared_axis_by_panel": axes,
        "sweep": {
            "n_pairings": int(len(sw)),
            "n_flagged_for_inspection": int(sw["flag_for_inspection"].sum()) if len(sw) else 0,
            "flagged": (sw.loc[sw["flag_for_inspection"], ["panel", "npi", "drug", "similarity",
                                                   "residual_z"]].to_dict("records")
                        if len(sw) else []),
            "flag": "top decile down both axes and residual z >= 3.0; an inspection "
                    "flag, not a significance threshold",
            "conclusion": None,
        },
        "npi_response_magnitude": {k: round(float(v), 4) for k, v in magnitude.items()},
        "rule": "A pairing is only specific if it ranks high down the drug axis AND "
                "down the NPI axis. A high rank on one axis alone is consistent with "
                "the drug merely sharing the panel's dominant axis with a low-magnitude "
                "NPI signature.",
    }
    if len(sw):
        verdict["sweep"]["calibration"] = CALIBRATION_NOTE
        print("\n" + CALIBRATION_NOTE)
        verdict["sweep"]["conclusion"] = (
            "No NPI-drug pairing is supported. The evidence is the two checks that need "
            "no null: the claimed pairing fails the reciprocal axis (erastin ranks "
            "cystine deprivation 3rd of 14 on the GEO consensus and 6th on the 51-cell-"
            "line Phase I consensus), and one axis carries 80% of the panel's variance, "
            "that axis being NPI response magnitude. The residual column below is "
            "descriptive only -- see the calibration note."
        )
        print("\n" + verdict["sweep"]["conclusion"])
    (out / "verdict.json").write_text(json.dumps(verdict, indent=2, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
