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

    verdict = {
        "pairs": rows, "shared_axis_by_panel": axes,
        "npi_response_magnitude": {k: round(float(v), 4) for k, v in magnitude.items()},
        "rule": "A pairing is only specific if it ranks high down the drug axis AND "
                "down the NPI axis. A high rank on one axis alone is consistent with "
                "the drug merely sharing the panel's dominant axis with a low-magnitude "
                "NPI signature.",
    }
    (out / "verdict.json").write_text(json.dumps(verdict, indent=2, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
