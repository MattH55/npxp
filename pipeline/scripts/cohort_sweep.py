"""Score every GSE95640 subject as a patient, holding each one out of the LCD signature.

For each of the 191 subjects:
  - the LCD signature is rebuilt from the other 190 pairs;
  - s_P is the subject's baseline z against the other 190 baselines, plus the
    obese-vs-lean offset from `npi-pharma build-offset`, i.e. a z against health;
  - every drug in drugs.parquet is scored with the LCD signature.

Reports how often LCD reverses s_P, which drugs rank first, whether the named
drug set outranks the random background (per-patient AUC), and how much the
per-patient rankings differ from one another. The shared offset is a large
common term, so strong agreement between patients is expected. That agreement
is reported, not hidden.

    python scripts/cohort_sweep.py --offset data/processed/offset_obese_vs_lean_GSE244118.tsv \
        --drugs data/processed/signatures/drugs.parquet --drug-set metabolic --out out/cohort_sweep
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from npi_pharma.config import load_config
from npi_pharma.gene_sets import load_gmt
from npi_pharma.ingest.catalog import _metadata, load_catalog
from npi_pharma.ingest.geo import counts_to_log_cpm, read_counts
from npi_pharma.interact.rank import rank_pairs
from npi_pharma.model import PROV_GEO
from npi_pharma.patient.encode import encode_patient
from npi_pharma.signatures.build import build_signature
from npi_pharma.store import load_signatures

NPI_ID = "LCD_adipose_GSE95640"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--raw-dir", default="data/raw")
    ap.add_argument("--catalog", default="configs/npi_catalog.yaml")
    ap.add_argument("--offset", required=True)
    ap.add_argument("--drugs", default="data/processed/signatures/drugs.parquet")
    ap.add_argument("--drug-sets", default="configs/drug_sets.yaml")
    ap.add_argument("--drug-set", default="metabolic")
    ap.add_argument("--limit", type=int, help="first N subjects only")
    ap.add_argument("--out", default="out/cohort_sweep")
    a = ap.parse_args(argv)

    cfg = load_config(None)
    gene_sets = load_gmt(None)
    entry = next(e for e in load_catalog(a.catalog, min_n=cfg["min_npi_n"]) if e.npi_id == NPI_ID)
    inp, raw = entry.record["inputs"], Path(a.raw_dir)
    idm = pd.read_csv(raw / inp["id_map"], sep="\t", header=None, index_col=0, dtype=str)[1]
    expr = counts_to_log_cpm(read_counts(raw / inp["counts"]), idm)
    sheet = pd.read_csv(raw / inp["samples"], sep="\t", index_col=0, dtype=str)
    pre_s = sheet[sheet["timepoint"] == "CID1"].reset_index().set_index("subject")["sample"]
    post_s = sheet[sheet["timepoint"] == "CID2"].reset_index().set_index("subject")["sample"]
    subjects = sorted(set(pre_s.index) & set(post_s.index))[: a.limit]
    pre = expr[pre_s.loc[sorted(pre_s.index)]].set_axis(sorted(pre_s.index), axis=1)
    post = expr[post_s.loc[sorted(post_s.index)]].set_axis(sorted(post_s.index), axis=1)
    baseline = expr[pre_s.loc[sorted(pre_s.index)].tolist()]

    off_df = pd.read_csv(a.offset, sep="\t", index_col=0)
    offset, off_name = off_df.iloc[:, 0], off_df.columns[0]
    drugs = load_signatures(a.drugs)
    named = set(yaml.safe_load(Path(a.drug_sets).read_text())["drug_sets"][a.drug_set])
    named = {"sirolimus" if d == "rapamycin" else d for d in named}
    md = _metadata(entry.record, PROV_GEO, entry.quality_flag)

    rows, ranks = [], {}
    for i, subj in enumerate(subjects):
        keep = [s for s in pre.columns if s != subj]
        npi = build_signature(pre[keep], post[keep], md | {"sample_size": len(keep)}, paired=True)
        patient = encode_patient(baseline, pre_s[subj], "adipose", offset=offset, offset_desc=off_name)
        table, _ = rank_pairs(patient, [npi], drugs, gene_sets, cfg)
        is_named = table["drug_id"].isin(named).to_numpy()
        c = table["composite"].to_numpy()
        x, y = c[is_named], c[~is_named]
        auc = float(np.mean([(u > v) + 0.5 * (u == v) for u in x for v in y])) if len(x) and len(y) else np.nan
        rows.append({
            "subject": subj, "sample": pre_s[subj], "reverse_npi": float(table["reverse_npi"].iloc[0]),
            "top_drug": table["drug_id"].iloc[0], "top_named_drug": table.loc[is_named, "drug_id"].iloc[0],
            "auc_named_vs_background": auc,
            "n_positive_complementarity": int((table["complementarity"] > 0).sum()),
        })
        ranks[subj] = table.set_index("drug_id")["composite_rank"]
        if (i + 1) % 25 == 0:
            print(f"  {i + 1}/{len(subjects)}", file=sys.stderr)

    res = pd.DataFrame(rows)
    rk = pd.DataFrame(ranks)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    res.to_csv(out / "per_patient.tsv", sep="\t", index=False)
    rk.to_csv(out / "drug_ranks.tsv", sep="\t")

    rho = rk.rank().corr(method="spearman").to_numpy()
    off_share = []
    for s in subjects[:25]:
        p = encode_patient(baseline, pre_s[s], "adipose", offset=offset, offset_desc=off_name)
        v = pd.Series(p.disease_vector, index=p.genes)
        o = offset.loc[v.index]
        off_share.append(float(np.dot(o, o) / np.dot(v, v)))

    print(f"subjects: {len(res)}  drugs: {len(drugs)} ({len(named & {d.sig_id for d in drugs})} in '{a.drug_set}')")
    print(f"LCD reverses s_P (reverse_npi > 0): {(res.reverse_npi > 0).mean():.1%}; "
          f"median reverse_npi {res.reverse_npi.median():+.3f}")
    print(f"AUC named set vs background: median {res.auc_named_vs_background.median():.3f}, "
          f"IQR {res.auc_named_vs_background.quantile(.25):.3f}-{res.auc_named_vs_background.quantile(.75):.3f}")
    print(f"patients with any positive complementarity: {(res.n_positive_complementarity > 0).mean():.1%}")
    print("top drug frequency: " + ", ".join(f"{k} {v}" for k, v in res.top_drug.value_counts().head(6).items()))
    print("top named drug frequency: " + ", ".join(f"{k} {v}" for k, v in res.top_named_drug.value_counts().head(6).items()))
    print(f"between-patient drug-rank Spearman: median {np.median(rho[np.triu_indices_from(rho, 1)]):.3f}")
    print(f"share of |s_P|^2 carried by the shared offset (first 25 patients): median {np.median(off_share):.2f}")
    print(f"mean composite rank of each drug:\n{rk.mean(axis=1).sort_values().round(1).to_string()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
