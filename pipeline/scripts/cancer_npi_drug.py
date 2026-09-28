"""NPI x drug sensitisation in cancer cell lines: build, gate, rank.

Predictor: an NPI sensitises cells to a drug when its expression signature moves
cells toward the baseline state of lines that are sensitive to that drug
(npi_pharma.cancer). Drug programs come from DepMap expression x PRISM 19Q4 AUC
(1,390 drugs), averaged with GDSC1/GDSC2 programs for drugs those screens share.
The primary score uses drug-specific programs, i.e. the general-sensitivity
axis (mean program over all drugs) projected out.

Pre-registered gates (fixed before any NPI was scored), from measured results:
  G1  heat applied to cancer cells in vitro scores platinum (cisplatin,
      oxaliplatin, carboplatin) above 5-FU and mitomycin C, which were additive
      with heat under the same protocol (Helderman 2020, tier-1 rows).
  G2  cystine deprivation x erastin scores positive and in the top 10% of drugs
      (tier-1: cystine withdrawal + erastin synergy).
  G3  methionine restriction x 5-FU and x oxaliplatin score positive (tier-2).
Not testable: glucose restriction x metformin (metformin is in neither PRISM
nor GDSC); stiffness (no signature); fasting in vivo (no cell-line signature).

Inputs:
  --cancer-resource  "Cancer Resource" folder of MattH55/Physiological-Fitness-Landscape (modifier_profiles)
  data/raw/depmap    Model.csv, OmicsExpressionProteinCodingGenesTPMLogp1.csv,
                     prism19q4_secondary_dose_response.csv, GDSC1/2 fitted dose response
  data/raw/GSE62673  npi-pharma fetch-geo GSE62673 (amino-acid deprivation, MCF7/PC3)
  data/raw/GSE70138  LINCS gene_info (landmark genes)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from npi_pharma.cancer.programs import (
    drug_specific, load_expression, load_gdsc_auc, load_prism_auc, resistance_programs, split_half_reliability,
)
from npi_pharma.cancer.sensitize import percentile, score_matrix
from npi_pharma.ingest.geo import collapse_probes, maybe_log2, read_series_matrix
from npi_pharma.signatures.build import build_signature

HEAT_IN_VITRO = [
    "gse48398_mcf7_heat45c30min", "gse48398_mda231_heat45c30min", "gse48398_mda468_heat45c30min",
    "gse48398_mcf10a_heat45c30min", "gse10043_u937_mildhyperthermia41c30min", "gse75127_hsc3_hyperthermia44c90min",
]
PLATINUM = ["cisplatin", "oxaliplatin", "carboplatin"]
HEAT_ADDITIVE = ["5-fluorouracil", "mitomycin-c"]
ANCHOR_DRUGS = PLATINUM + HEAT_ADDITIVE + ["erastin", "doxorubicin", "paclitaxel", "temozolomide", "cyclophosphamide"]
IN_VITRO_CONTEXTS = ("mcf7", "t47d", "mdamb231", "mdamb468", "mcf10a", "u937", "hsc3", "u87", "lovo")


def zcol(df: pd.DataFrame) -> pd.DataFrame:
    return df.apply(lambda c: (c - c.mean()) / c.std(ddof=0))


def gse62673_signatures(raw: Path) -> pd.DataFrame:
    """Cystine and methionine deprivation vs control media, 24 h (MCF7; PC3 methionine)."""
    expr, s = read_series_matrix(raw / "GSE62673" / "GSE62673_series_matrix.txt.gz")
    pm = pd.read_csv(raw / "GSE62673" / "probe_map.tsv", sep="\t", header=None, index_col=0, dtype=str)[1]
    expr = maybe_log2(collapse_probes(expr, pm))
    s = s[s["time"] == "24 hours"]
    out = {}
    arms = {
        "gse62673_mcf7_cystine_deprivation24h": ("MCF7", "deprived of Cystine", "control media"),
        "gse62673_mcf7_methionine_deprivation24h": ("MCF7", "deprived of Methionine", "control media"),
        "gse62673_pc3_methionine_0um24h": ("PC3", "0 μM Methionine media", "20 μM Methionine media"),
    }
    for sid, (cell, trt, ctrl) in arms.items():
        post = s.index[(s["cell line"] == cell) & (s["treatment"] == trt)]
        pre = s.index[(s["cell line"] == cell) & (s["treatment"] == ctrl)]
        if len(post) < 2 or len(pre) < 2:
            raise ValueError(f"{sid}: {len(pre)} control / {len(post)} treated samples")
        sig = build_signature(expr[pre], expr[post], {"npi_id": sid, "contrast": f"{trt} vs {ctrl}"}, paired=False)
        out[sid] = sig.as_series()
    return pd.DataFrame(out)


def modifier_profiles(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    sigs, meta = {}, {}
    for d in sorted((root / "data" / "modifier_profiles").iterdir()):
        m = pd.read_csv(d / "migep.csv", sep=None, engine="python", index_col=0).iloc[:, 0]
        m.index = m.index.str.upper()
        sigs[d.name] = m[~m.index.duplicated()]
        md = json.loads((d / "metadata.json").read_text())
        meta[d.name] = {"context": md.get("biological_context"), "tissue": md.get("tissue"),
                        "in_vitro_cancer_cell": str(md.get("biological_context", "")).lower() in IN_VITRO_CONTEXTS}
    return pd.DataFrame(sigs), pd.DataFrame(meta).T


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cancer-resource", required=True)
    ap.add_argument("--raw-dir", default="data/raw")
    ap.add_argument("--out", default="out/cancer_npi_drug")
    ap.add_argument("--top", type=int, default=15)
    ap.add_argument("--rebuild", action="store_true", help="ignore the cached drug programs")
    a = ap.parse_args(argv)
    raw, dm, out = Path(a.raw_dir), Path(a.raw_dir) / "depmap", Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    cache = out / "programs_cache.pkl"
    if cache.exists() and not a.rebuild:
        programs, specific, sources, n_prism, rel = pd.read_pickle(cache)
    else:
        gi = pd.read_csv(next((raw / "GSE70138").glob("*gene_info*.txt")), sep="\t", dtype=str)
        landmarks = sorted(set(gi.loc[gi["pr_is_lm"] == "1", "pr_gene_symbol"].str.upper()))
        models = pd.read_csv(dm / "Model.csv")
        lineage = models.set_index("ModelID")["OncotreeLineage"]
        expr = load_expression(dm / "OmicsExpressionProteinCodingGenesTPMLogp1.csv", landmarks + ["SLFN11", "ABCB1"])
        print(f"expression: {expr.shape[0]} lines x {expr.shape[1]} genes", file=sys.stderr)

        prism = load_prism_auc(dm / "prism19q4_secondary_dose_response.csv")
        p_prism, n_prism = resistance_programs(expr, prism, lineage)
        parts = {"PRISM": p_prism}
        for name in ("GDSC1", "GDSC2"):
            g = load_gdsc_auc(dm / f"{name}_fitted_dose_response_25Feb20.csv", models)
            parts[name], _ = resistance_programs(expr, g[[c for c in g if c in ANCHOR_DRUGS]], lineage)
        programs = zcol(p_prism.copy())
        sources = {d: ["PRISM"] for d in programs}
        for d in ANCHOR_DRUGS:
            avail = [zcol(parts[k][[d]])[d] for k in parts if d in parts[k]]
            if avail:
                programs[d] = pd.concat(avail, axis=1).mean(axis=1)
                sources[d] = [k for k in parts if d in parts[k]]
        programs = programs.drop(index=["SLFN11", "ABCB1"], errors="ignore")
        specific = drug_specific(programs)

        rel = {d: split_half_reliability(expr, prism[d], lineage, n_splits=6) for d in ANCHOR_DRUGS if d in prism}
        pd.to_pickle((programs, specific, sources, n_prism, rel), cache)

    mods, meta = modifier_profiles(Path(a.cancer_resource))
    aa = gse62673_signatures(raw)
    npis = pd.concat([mods, aa.reindex(mods.index.union(aa.index))], axis=1)
    for c in aa:
        meta.loc[c] = {"context": c.split("_")[1], "tissue": c.split("_")[1].upper(), "in_vitro_cancer_cell": True}

    scores = {"specific": score_matrix(npis, specific), "raw": score_matrix(npis, programs)}
    pct = {k: percentile(v) for k, v in scores.items()}

    gates = {}
    for k in ("specific", "raw"):
        s, p = scores[k], pct[k]
        heat = s.loc[HEAT_IN_VITRO]
        plat = heat[[d for d in PLATINUM if d in heat]].mean(axis=1)
        addv = heat[[d for d in HEAT_ADDITIVE if d in heat]].mean(axis=1)
        cys = "gse62673_mcf7_cystine_deprivation24h"
        mets = ["gse62673_mcf7_methionine_deprivation24h", "gse62673_pc3_methionine_0um24h"]
        g = {
            "G1_heat_platinum_minus_additive": {
                "per_signature": (plat - addv).round(4).to_dict(),
                "mean_difference": float((plat - addv).mean()),
                "n_signatures_platinum_higher": int(((plat - addv) > 0).sum()),
                "pass": bool((plat - addv).mean() > 0 and ((plat - addv) > 0).sum() >= 4),
            },
            "G2_cystine_erastin": {
                "score": float(s.loc[cys, "erastin"]), "percentile": float(p.loc[cys, "erastin"]),
                "pass": bool(s.loc[cys, "erastin"] > 0 and p.loc[cys, "erastin"] >= 0.9),
            },
            "G3_methionine_5fu_oxaliplatin": {
                "scores": {m: {d: float(s.loc[m, d]) for d in ("5-fluorouracil", "oxaliplatin")} for m in mets},
                "pass": bool(all(s.loc[m, d] > 0 for m in mets for d in ("5-fluorouracil", "oxaliplatin"))),
            },
        }
        g["all_pass"] = all(v["pass"] for v in g.values() if isinstance(v, dict) and "pass" in v)
        gates[k] = g

    s = scores["specific"]
    long = s.stack().rename("score").reset_index().rename(columns={"level_0": "npi", "level_1": "drug"})
    long.columns = ["npi", "drug", "score"]
    long["percentile"] = pct["specific"].stack().to_numpy()
    long["in_vitro_cancer_cell_npi"] = long["npi"].map(meta["in_vitro_cancer_cell"]).fillna(False).astype(bool)
    long["drug_program_sources"] = long["drug"].map(lambda d: "+".join(sources.get(d, ["PRISM"])))
    long["drug_n_lines_prism"] = long["drug"].map(n_prism)
    long = long.sort_values("score", ascending=False)
    long.to_csv(out / "npi_drug_scores.tsv", sep="\t", index=False)

    summary = {
        "n_npis": int(npis.shape[1]), "n_drugs": int(s.shape[1]), "n_genes": int(len(npis.index.intersection(specific.index))),
        "gates": gates, "prism_split_half_reliability": {k: round(v, 3) for k, v in rel.items()},
        "anchor_program_sources": {d: sources[d] for d in ANCHOR_DRUGS if d in sources},
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    for k in ("specific", "raw"):
        g = gates[k]
        print(f"[{k}] G1 heat platinum-minus-additive mean {g['G1_heat_platinum_minus_additive']['mean_difference']:+.4f} "
              f"({g['G1_heat_platinum_minus_additive']['n_signatures_platinum_higher']}/6 higher) "
              f"{'PASS' if g['G1_heat_platinum_minus_additive']['pass'] else 'FAIL'}; "
              f"G2 cystine x erastin {g['G2_cystine_erastin']['score']:+.4f} pct {g['G2_cystine_erastin']['percentile']:.3f} "
              f"{'PASS' if g['G2_cystine_erastin']['pass'] else 'FAIL'}; "
              f"G3 methionine {g['G3_methionine_5fu_oxaliplatin']['scores']} "
              f"{'PASS' if g['G3_methionine_5fu_oxaliplatin']['pass'] else 'FAIL'}; all: {g['all_pass']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
