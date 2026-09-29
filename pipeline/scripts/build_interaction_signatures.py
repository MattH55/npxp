"""Build measured interaction signatures from the factorial GEO series found.

Reads `out/combination_series/factorial_series.tsv` (from
scripts/find_combination_series.py), loads each series' expression, and computes

    I(A, B) = combo - A - B + control

the transcriptional excess over additivity (npi_pharma.interact.measured). This is
observed interaction data, which is what every inferred approach in this project
was missing: signature composition carries no detectable information about synergy
across 50k drug-pair observations (docs/validation_drugcomb.md).

A series is used only when all four arms resolve to samples and the expression
columns map to those samples through GEO's own labels. Nothing is guessed.

    python scripts/build_interaction_signatures.py
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from build_drug_consensus import load_symbol_reference, supplementary_table  # noqa: E402
from find_combination_series import detect_factorial  # noqa: E402

from npi_pharma.ingest.fetch import fetch_geo_series  # noqa: E402
from npi_pharma.ingest.geo import collapse_probes, maybe_log2, read_series_matrix  # noqa: E402
from npi_pharma.interact.measured import (  # noqa: E402
    additivity_summary, build_interaction_signature, interaction_contrast,
)
from npi_pharma.store import save_signatures  # noqa: E402


def load_expression(gse: str, raw: Path, samples: pd.DataFrame, need: list[str],
                    symbols: set[str] | None) -> pd.DataFrame | None:
    """Genes x GSM on a log scale, from the series matrix or a supplementary table."""
    mats = sorted((raw / gse).glob("*series_matrix.txt.gz"))
    expr, _ = read_series_matrix(mats[0]) if mats else (pd.DataFrame(), None)
    if expr.shape[0] >= 500:
        pm = raw / gse / "probe_map.tsv"
        if not pm.exists():
            alt = sorted((raw / gse).glob("probe_map_*.tsv"))
            pm = alt[0] if alt else None
        if pm is not None and pm.exists():
            m = pd.read_csv(pm, sep="\t", header=None, index_col=0, dtype=str)[1]
            expr = collapse_probes(expr, m)
        else:
            expr.index = expr.index.astype(str).str.upper()
            expr = expr.groupby(level=0).mean()
        return maybe_log2(expr)
    counts = supplementary_table(raw / gse, samples, need=need, symbols=symbols, why=[])
    if counts is None:
        return None
    total = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts / total * 1e6
    keep = (cpm >= 1).mean(axis=1) >= 0.5
    return np.log2(cpm[keep] + 1)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--factorial", default="out/combination_series/factorial_series.tsv")
    ap.add_argument("--raw-dir", default="data/raw")
    ap.add_argument("--min-genes", type=int, default=3000)
    ap.add_argument("--out", default="data/processed/signatures/interactions_measured.parquet")
    ap.add_argument("--report", default="out/interaction_signatures")
    a = ap.parse_args(argv)
    raw, rep = Path(a.raw_dir), Path(a.report)
    rep.mkdir(parents=True, exist_ok=True)
    symbols = load_symbol_reference(raw)

    table = pd.read_csv(a.factorial, sep="\t")
    sigs, log = [], []
    for r in table.itertuples():
        gse = r.accession
        try:
            fetch_geo_series(gse, raw, suppl=True, log=lambda *x, **k: None)
            mats = sorted((raw / gse).glob("*series_matrix.txt.gz"))
            _, samples = read_series_matrix(mats[0])
            hit = detect_factorial(samples)
            if hit is None or len(hit["singles"]) < 2:
                log.append({"accession": gse, "outcome": "factorial design no longer resolves"})
                continue
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
            if expr is None:
                log.append({"accession": gse, "outcome": "no usable expression table"})
                continue
            missing = [s for s in need if s not in expr.columns]
            if missing:
                log.append({"accession": gse, "outcome": f"{len(missing)} arm samples absent from the table"})
                continue
            if expr.shape[0] < a.min_genes:
                log.append({"accession": gse, "outcome": f"only {expr.shape[0]} genes"})
                continue
            meta = {"sig_id": f"{a_name} + {b_name}|{gse}", "agent_a": a_name, "agent_b": b_name,
                    "cell_line": None, "source_accessions": [gse], "field": hit["field"],
                    "title": getattr(r, "title", "")}
            sig = build_interaction_signature(expr[arms["control"]], expr[arms["a"]],
                                              expr[arms["b"]], expr[arms["combo"]], meta)
            effect, _ = interaction_contrast(expr[arms["control"]], expr[arms["a"]],
                                             expr[arms["b"]], expr[arms["combo"]])
            summ = additivity_summary(effect)
            sigs.append(sig)
            log.append({"accession": gse, "outcome": "built", "agent_a": a_name, "agent_b": b_name,
                        "n_genes": len(sig.genes), "quality_flag": sig.quality_flag,
                        "sd_of_interaction": round(summ["sd_of_interaction"], 3),
                        "n_per_arm": json.dumps({k: len(v) for k, v in arms.items()})})
            print(f"  {gse}: {a_name} + {b_name} -> {len(sig.genes)} genes, "
                  f"sd {summ['sd_of_interaction']:.3f}", file=sys.stderr)
        except Exception as e:
            log.append({"accession": gse, "outcome": f"error: {type(e).__name__}: {e}"})
    pd.DataFrame(log).to_csv(rep / "log.tsv", sep="\t", index=False)
    if sigs:
        save_signatures(sigs, a.out)
    built = [x for x in log if x.get("outcome") == "built"]
    print(f"\n{len(built)} measured interaction signatures from {len(table)} factorial series "
          f"-> {a.out}")
    if built:
        print(pd.DataFrame(built)[["accession", "agent_a", "agent_b", "n_genes",
                                   "sd_of_interaction", "quality_flag"]].to_string(index=False))
    fails = pd.DataFrame([x for x in log if x.get("outcome") != "built"])
    if len(fails):
        print("\nnot built:")
        print(fails["outcome"].str.slice(0, 50).value_counts().to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
