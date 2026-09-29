"""Single-agent drug signatures from the factorial series' own control arms.

The factorial corpus was collected for the interaction contrast, which turned out
not to reproduce across independent series (docs/measured_interactions.md). The same
four arms also give each agent's *main effect*, and those do reproduce: in the one
replicate pair available, abemaciclib's main effect agrees at cosine 0.512 and the
combination's at 0.632, against 0.18 for the same drug in two cell lines elsewhere
in this project (docs/npi_drug_retrieval.md). So the main effects are worth keeping
even though the interaction from the identical data is not.

Each series contributes two signatures, ``a - control`` and ``b - control``, built
the same way as every other drug signature in the panel (moderated statistic, then
standardised) so they are directly comparable. Two of them are also a free check on
the drug search: GSE338229's cisplatin arm is a drug the search found independently,
so the two routes can be compared.

These are single-series, single-cell-line signatures, which this project measured to
be dominated by the cell line rather than the drug. They are therefore labelled
"very low" reliability by the same rule as every other single-series signature
(scripts/rank_npi_drug_all.py) and are useful as consensus *inputs*, not on their own.

    python scripts/drug_signatures_from_factorial.py
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from build_drug_consensus import load_symbol_reference  # noqa: E402
from build_interaction_signatures import load_expression  # noqa: E402
from find_combination_series import detect_factorial  # noqa: E402

from npi_pharma.ingest.geo import read_series_matrix  # noqa: E402
from npi_pharma.model import DRUG  # noqa: E402
from npi_pharma.signatures.build import build_signature  # noqa: E402
from npi_pharma.store import save_signatures  # noqa: E402


# A drug-resistant derivative line is named by a suffix: IGROV-1/CP is IGROV-1 made
# cisplatin-resistant, and the response of a resistant line to the drug it resists is
# not that drug's response. build_drug_consensus.py rejects resistant *arms*
# (RESISTANCE_VALUE) but says nothing about the cell line, which is annotated
# elsewhere in the series. Flagged here rather than dropped: the signature is real,
# it just answers a different question.
RESISTANT_LINE = re.compile(
    r"(?:^|[/\-_ ])(?:CP|CDDP|DDP|CisR|CIS|R|RES|MDR|TXR|DOXR|PTXR)(?:$|[/\-_ 0-9])|"
    r"resistant|refractory", re.I)


def resistant_line_flag(cell_line: str | None) -> str | None:
    """"resistant_cell_line" when the line's name says it was selected for resistance."""
    if not cell_line:
        return None
    return "resistant_cell_line" if RESISTANT_LINE.search(str(cell_line)) else None


def agent_label(raw: str) -> str:
    """The drug name out of a GEO arm label, or the label itself if nothing is clearer.

    GEO arm values carry dose and time ("Dalpiciclib 2uM, 48h") and parenthetical
    synonyms ("ABT199 (venetoclax)"). Only the obvious decorations are stripped --
    the raw label is kept in the signature's metadata either way, so nothing is lost
    if this guesses wrong.
    """
    s = str(raw).split(",")[0].strip()
    s = s.split("(")[0].strip() if "(" in s else s
    for cut in (" treatment at ", " treatment", " at "):
        if cut in s:
            s = s.split(cut)[0].strip()
    # a trailing dose: "Enzalutamide 35uM", "A51 62.5nM"
    parts = s.split()
    while len(parts) > 1 and any(u in parts[-1].lower() for u in ("um", "nm", "mm", "μm", "µm", "h")):
        parts = parts[:-1]
    return " ".join(parts) or str(raw).strip()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--factorial", default="out/combination_series/factorial_series.tsv")
    ap.add_argument("--raw-dir", default="data/raw")
    ap.add_argument("--min-genes", type=int, default=3000)
    ap.add_argument("--out", default="data/processed/signatures/drugs_from_factorial.parquet")
    ap.add_argument("--report", default="out/drug_signatures_factorial")
    a = ap.parse_args(argv)
    raw, rep = Path(a.raw_dir), Path(a.report)
    rep.mkdir(parents=True, exist_ok=True)
    symbols = load_symbol_reference(raw)

    sigs, log = [], []
    for r in pd.read_csv(a.factorial, sep="\t").itertuples():
        gse = r.accession
        try:
            mats = sorted((raw / gse).glob("*series_matrix.txt.gz"))
            if not mats:
                log.append({"accession": gse, "outcome": "series matrix not downloaded"})
                continue
            _, samples = read_series_matrix(mats[0])
            hit = detect_factorial(samples)
            if hit is None or len(hit["singles"]) < 2:
                log.append({"accession": gse, "outcome": "factorial design no longer resolves"})
                continue
            vals = samples[hit["field"]].astype(str)
            control = vals.index[vals.isin(hit["control"])].tolist()
            arms = {s: vals.index[vals == s].tolist() for s in hit["singles"][:2]}
            need = control + [s for v in arms.values() for s in v]
            expr = load_expression(gse, raw, samples, need, symbols)
            if expr is None or any(s not in expr.columns for s in need):
                log.append({"accession": gse, "outcome": "no usable expression table"})
                continue
            if expr.shape[0] < a.min_genes:
                log.append({"accession": gse, "outcome": f"only {expr.shape[0]} genes"})
                continue
            cell = samples["source_name_ch1"].iloc[0] if "source_name_ch1" in samples else None
            for label, gsms in arms.items():
                drug = agent_label(label)
                # Built exactly as scripts/build_drug_consensus.py builds a panel
                # signature, so the two are directly comparable and can feed one
                # consensus: build_signature is the NPI constructor, and the drug
                # panel re-labels its output rather than duplicating the statistics.
                sig = build_signature(
                    expr[control], expr[gsms],
                    {"npi_id": f"{drug}|{gse}", "source_accessions": [gse],
                     "provenance": "geo_factorial_main_effect", "tissue": cell,
                     "contrast": f"{label} vs control ({hit['field']}), unpaired; "
                                 "single agent of a factorial design"},
                    paired=False, min_n=2)
                sig.kind = DRUG
                sig.modality = "compound"
                res = resistant_line_flag(cell)
                sig.meta = dict(sig.meta or {}) | {
                    "drug_id": drug, "arm_label": label, "cell_line": cell,
                    "field": hit["field"], "n_control": len(control), "n_treated": len(gsms),
                    "cell_line_flag": res}
                if res:
                    sig.quality_flag = f"{sig.quality_flag},{res}" if sig.quality_flag else res
                sigs.append(sig)
                log.append({"accession": gse, "outcome": "built", "drug": drug,
                            "arm_label": label, "cell_line": cell, "n_genes": len(sig.genes),
                            "n_treated": len(gsms), "n_control": len(control),
                            "quality_flag": sig.quality_flag, "cell_line_flag": res})
                print(f"  {gse}: {drug} ({len(gsms)}v{len(control)}) -> {len(sig.genes)} genes",
                      file=sys.stderr)
        except Exception as e:
            log.append({"accession": gse, "outcome": f"error: {type(e).__name__}: {e}"})
    t = pd.DataFrame(log)
    t.to_csv(rep / "log.tsv", sep="\t", index=False)
    if sigs:
        save_signatures(sigs, a.out)
    built = t[t["outcome"] == "built"] if len(t) else t
    (rep / "summary.json").write_text(json.dumps(
        {"n_signatures": len(sigs), "n_series": int(built["accession"].nunique()) if len(built) else 0,
         "drugs": sorted(built["drug"].unique()) if len(built) else [],
         "reliability": "very low -- single series, single cell line; useful as consensus "
                        "inputs, not on their own"}, indent=2))
    print(f"\n{len(sigs)} single-agent signatures from "
          f"{built['accession'].nunique() if len(built) else 0} factorial series -> {a.out}")
    if len(built):
        print(built[["accession", "drug", "cell_line", "n_genes", "n_treated",
                     "quality_flag"]].to_string(index=False))
    fails = t[t["outcome"] != "built"] if len(t) else t
    if len(fails):
        print("\nnot built:")
        print(fails["outcome"].str.slice(0, 48).value_counts().to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
