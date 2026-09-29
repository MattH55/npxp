"""Does a MEASURED interaction signature reproduce across independent series?

The measured interaction contrast (npi_pharma.interact.measured) was built because
every *inferred* interaction score in this project failed: signature composition
carries no detectable information about synergy across 50k drug-pair observations
(docs/validation_drugcomb.md). Computing the interaction directly from a factorial
design is different in kind -- it is observed, not predicted. But "observed" is not
the same as "reproducible", and that has to be measured too rather than assumed.

The corpus contains two independent series that profile the SAME combination:
abemaciclib + fulvestrant in MCF7 (GSE336734) and in CAMA1 (GSE336729). That makes
a direct test possible, and it comes with its own control. From the same four arms
of the same two series, this compares:

  * each single agent's main effect      (drug - control)
  * the combination's main effect        (combo - control)
  * the interaction contrast             (combo - a - b + control)

If the main effects reproduce and the interaction does not, the failure is specific
to the interaction term -- four group means differenced against each other -- and
not to the data, the download path or the gene mapping, all of which are shared.
That is the point of running the main effects as a positive control: it makes the
bar a fair one.

Cross-cell-line agreement is itself limited: the same drug in two cell lines agrees
at only cosine 0.18 in this project's own measurements (docs/npi_drug_retrieval.md),
so the main effects set the ceiling the interaction should be judged against, not
1.0.

    python scripts/validate_interaction_reproducibility.py
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from build_drug_consensus import load_symbol_reference  # noqa: E402
from build_interaction_signatures import load_expression  # noqa: E402
from find_combination_series import detect_factorial  # noqa: E402

from npi_pharma.ingest.fetch import fetch_geo_series  # noqa: E402
from npi_pharma.ingest.geo import read_series_matrix  # noqa: E402
from npi_pharma.interact.measured import interaction_contrast  # noqa: E402

# Series pairs that profile the same combination independently. Add a pair here only
# when both series really are separate submissions; two arms of one accession share
# protocol, batch and operator and would measure repeatability, not reproducibility.
REPLICATE_PAIRS = [("GSE336734", "GSE336729", "abemaciclib + fulvestrant")]

CONTRASTS = ["agent_a_main", "agent_b_main", "combination_main", "interaction"]


def cosine(x: pd.Series, y: pd.Series) -> tuple[float, int]:
    """Cosine between two contrasts on their shared genes, each standardised."""
    g = x.dropna().index.intersection(y.dropna().index)
    if len(g) < 200:
        return float("nan"), len(g)
    a, b = x[g], y[g]
    za = (a - a.mean()) / a.std(ddof=0)
    zb = (b - b.mean()) / b.std(ddof=0)
    return float((za * zb).mean()), len(g)


def contrasts_for(gse: str, raw: Path, symbols: set[str] | None) -> dict | None:
    """The four contrasts above from one factorial series, plus what it profiled."""
    fetch_geo_series(gse, raw, suppl=True, log=lambda *x, **k: None)
    mats = sorted((raw / gse).glob("*series_matrix.txt.gz"))
    if not mats:
        return None
    _, samples = read_series_matrix(mats[0])
    hit = detect_factorial(samples)
    if hit is None or len(hit["singles"]) < 2:
        return None
    vals = samples[hit["field"]].astype(str)
    a_name, b_name = hit["singles"][0], hit["singles"][1]
    arms = {
        "control": vals.index[vals.isin(hit["control"])].tolist(),
        "a": vals.index[vals == a_name].tolist(),
        "b": vals.index[vals == b_name].tolist(),
        "combo": vals.index[vals.isin(hit["combination"])].tolist(),
    }
    need = [s for v in arms.values() for s in v]
    expr = load_expression(gse, raw, samples, need, symbols)
    if expr is None or any(s not in expr.columns for s in need):
        return None

    def mean(key: str) -> pd.Series:
        return expr[arms[key]].mean(axis=1)

    effect, _ = interaction_contrast(expr[arms["control"]], expr[arms["a"]],
                                     expr[arms["b"]], expr[arms["combo"]])
    src = samples["source_name_ch1"].iloc[0] if "source_name_ch1" in samples else None
    return {
        "accession": gse, "agent_a": a_name, "agent_b": b_name, "cell_line": src,
        "n_per_arm": {k: len(v) for k, v in arms.items()},
        "agent_a_main": mean("a") - mean("control"),
        "agent_b_main": mean("b") - mean("control"),
        "combination_main": mean("combo") - mean("control"),
        "interaction": effect,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--raw-dir", default="data/raw")
    ap.add_argument("--out", default="out/interaction_reproducibility")
    a = ap.parse_args(argv)
    raw, out = Path(a.raw_dir), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    symbols = load_symbol_reference(raw)

    rows, notes = [], []
    for left, right, label in REPLICATE_PAIRS:
        L, R = contrasts_for(left, raw, symbols), contrasts_for(right, raw, symbols)
        if L is None or R is None:
            notes.append(f"{label}: {left if L is None else right} did not resolve")
            continue
        print(f"== {label}")
        for side in (L, R):
            print(f"   {side['accession']}  {side['cell_line']}  "
                  f"{side['agent_a']} / {side['agent_b']}  n={side['n_per_arm']}")
        for c in CONTRASTS:
            cos, n = cosine(L[c], R[c])
            rows.append({"pair": label, "left": left, "right": right, "contrast": c,
                         "cosine": round(cos, 4), "n_genes": n,
                         "left_cell_line": L["cell_line"], "right_cell_line": R["cell_line"]})
            print(f"   {c:18s} cos {cos:+.3f}  over {n} genes")

    t = pd.DataFrame(rows)
    t.to_csv(out / "reproducibility.tsv", sep="\t", index=False)
    verdict = None
    if len(t):
        main_eff = t[t["contrast"] != "interaction"]["cosine"]
        inter = t[t["contrast"] == "interaction"]["cosine"]
        verdict = {
            "n_pairs": int(t["pair"].nunique()),
            "main_effect_cosine_median": float(main_eff.median()),
            "main_effect_cosine_min": float(main_eff.min()),
            "interaction_cosine_median": float(inter.median()),
            "reproducible": bool(inter.median() >= main_eff.min()),
            "notes": notes,
            "caveat": "A single independent pair, at n=2 and n=3 per arm, in two "
                      "different cell lines. The main effects are the positive "
                      "control: they share every step of the pipeline with the "
                      "interaction, so a split between them is specific to the "
                      "interaction contrast. One pair cannot say how often that "
                      "happens, only that it does here.",
        }
        print(f"\nmain effects reproduce at median cos {verdict['main_effect_cosine_median']:+.3f} "
              f"(worst {verdict['main_effect_cosine_min']:+.3f}); "
              f"interaction at {verdict['interaction_cosine_median']:+.3f}")
        print("VERDICT: the measured interaction contrast does "
              f"{'' if verdict['reproducible'] else 'NOT '}reproduce across "
              "independent series at these sample sizes.")
    (out / "verdict.json").write_text(json.dumps(verdict, indent=2))
    for n in notes:
        print(f"  note: {n}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
