"""Multi-cell-line consensus drug signatures from several independent GEO series.

Why: a single-cell-line drug signature is dominated by the cell line rather than
the drug. Measured in this project (docs/npi_drug_retrieval.md): two *different*
DNA-damaging drugs in one cell line agree at cosine 0.71, while the *same* drug in
two cell lines agrees at only 0.18. LINCS avoids this by taking a median over many
cell lines, which is why the one validated retrieval result used LINCS consensus.
LINCS Phase II, though, lacks every platinum agent, 5-FU and erastin, and Phase I
(21 GB) does not fit here. So a consensus has to be assembled from independent GEO
series.

Method: search GEO per drug, keep series whose *series matrix carries the
expression table* (so sample columns are GSMs that GEO's own annotation labels,
and no column name has to be guessed), auto-detect a control arm and a
drug-treated arm from the sample annotations, build one signature per series, then
take the per-gene median across series as the consensus.

Validation built in: for drugs present in BOTH LINCS and GEO (doxorubicin,
paclitaxel, mitomycin C, temozolomide, metformin) the GEO consensus is compared
with the LINCS consensus. If the method recovers LINCS where both exist, the
consensus it builds for the platinums LINCS lacks is trustworthy on the same terms.

    python scripts/build_drug_consensus.py --drugs cisplatin oxaliplatin carboplatin
    python scripts/build_drug_consensus.py --validate-against data/processed/signatures/lincs_all.parquet
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from npi_pharma.ingest.fetch import fetch_geo_series
from npi_pharma.ingest.geo import collapse_probes, maybe_log2, read_series_matrix
from npi_pharma.model import DRUG
from npi_pharma.signatures.build import build_signature, standardize
from npi_pharma.store import load_signatures, save_signatures

sys.path.insert(0, str(Path(__file__).parent))
from find_drug_datasets import esearch, esummary, score  # noqa: E402

CONTROL = re.compile(r"\b(control|untreated|vehicle|dmso|mock|no treatment|none|normal|"
                     r"non.?treated|nt|ctrl|wt|parental|0 ?(u|µ|μ)m|baseline)\b", re.I)
# a treated arm must name the drug; aliases keep common spellings together
ALIASES = {
    "cisplatin": ["cisplatin", "cddp", "cis-platin", "cis platinum"],
    "carboplatin": ["carboplatin", "cbdca"],
    "oxaliplatin": ["oxaliplatin", "l-ohp", "oxaliplatinum"],
    "5-fluorouracil": ["5-fluorouracil", "5-fu", "fluorouracil", "5fu"],
    "mitomycin C": ["mitomycin"],
    "doxorubicin": ["doxorubicin", "adriamycin", "dox"],
    "paclitaxel": ["paclitaxel", "taxol"],
    "temozolomide": ["temozolomide", "tmz"],
    "metformin": ["metformin"],
    "erastin": ["erastin"],
    "triapine": ["triapine", "3-ap"],
    "lomustine": ["lomustine", "ccnu"],
    "cyclophosphamide": ["cyclophosphamide", "ctx", "4-hydroxycyclophosphamide", "mafosfamide"],
    # Further widely used drugs that LINCS Phase II also lacks entirely.
    "methotrexate": ["methotrexate", "mtx", "amethopterin"],
    "pemetrexed": ["pemetrexed", "alimta"],
    "topotecan": ["topotecan", "hycamtin"],
    "idarubicin": ["idarubicin"],
    "vincristine": ["vincristine", "oncovin"],
    "vinblastine": ["vinblastine"],
    "vinorelbine": ["vinorelbine", "navelbine"],
    "melphalan": ["melphalan", "l-pam"],
    "carmustine": ["carmustine", "bcnu"],
    "bleomycin": ["bleomycin"],
    "actinomycin D": ["actinomycin d", "actinomycin-d", "dactinomycin"],
    "5-azacytidine": ["5-azacytidine", "azacitidine", "5-aza", "5-azacitidine"],
    # Present in LINCS: used to validate the method, not because a consensus is needed.
    "etoposide": ["etoposide", "vp-16", "vp16"],
    "gemcitabine": ["gemcitabine"],
    "vorinostat": ["vorinostat", "saha", "suberoylanilide"],
    "tamoxifen": ["tamoxifen", "4-hydroxytamoxifen", "4-oht"],
    "bortezomib": ["bortezomib", "ps-341", "velcade"],
}
# In LINCS Phase II, so a GEO consensus can be checked against a trusted one.
IN_LINCS_FOR_VALIDATION = ["doxorubicin", "paclitaxel", "temozolomide", "metformin", "mitomycin C",
                           "etoposide", "gemcitabine", "vorinostat", "tamoxifen", "bortezomib"]
# Drugs LINCS Phase II does not carry, so a GEO consensus is the only route to a
# multi-cell-line signature for them. `cetuximab` is excluded (an antibody, not in
# compound screens) and `rapamycin` is covered by LINCS under the name `sirolimus`.
MISSING_FROM_LINCS = ["cisplatin", "carboplatin", "oxaliplatin", "5-fluorouracil", "erastin",
                      "triapine", "lomustine", "cyclophosphamide", "methotrexate", "pemetrexed",
                      "topotecan", "idarubicin", "vincristine", "vinblastine", "vinorelbine",
                      "melphalan", "carmustine", "bleomycin", "actinomycin D", "5-azacytidine"]
# A resistant-vs-parental comparison measures adaptation, not drug response, and its
# annotations look deceptively like a control/treated split ("parental" vs "resistant").
RESISTANCE_FIELD = re.compile(r"resistan|sensitiv|subline|clone", re.I)
RESISTANCE_VALUE = re.compile(r"resistant|refractory|\bR\d*\b$|parental", re.I)
SKIP_FIELD = {"title", "geo_accession", "status", "submission_date", "last_update_date", "type",
              "channel_count", "molecule_ch1", "extract_protocol_ch1", "taxid_ch1", "data_processing",
              "platform_id", "contact_name", "contact_email", "contact_institute", "contact_address",
              "contact_city", "contact_country", "supplementary_file", "supplementary_file_1",
              "relation", "data_row_count", "instrument_model", "library_selection", "library_source",
              "library_strategy", "description", "organism_ch1", "label_ch1", "hyb_protocol",
              "scan_protocol", "contact_laboratory", "contact_department", "contact_zip/postal_code"}
# combination arms must not be mistaken for the drug alone
COMBO = re.compile(r"\+|plus |combination|combined|co.?treat|with ", re.I)


def detect_arms(samples: pd.DataFrame, drug: str) -> tuple[list[str], list[str], str] | None:
    """Find (control GSMs, treated GSMs, field) from GEO's sample annotations."""
    names = ALIASES.get(drug, [drug])
    pat = re.compile("|".join(re.escape(n) for n in names), re.I)
    best = None
    for field in samples.columns:
        if field in SKIP_FIELD or RESISTANCE_FIELD.search(field):
            continue
        vals = samples[field].astype(str)
        if vals.nunique() < 2 or vals.nunique() > 12:
            continue
        if any(RESISTANCE_VALUE.search(v) for v in vals.unique()):
            continue                      # adaptation comparison, not acute treatment
        ctrl_vals = [v for v in vals.unique() if CONTROL.search(v) and not pat.search(v)]
        drug_vals = [v for v in vals.unique() if pat.search(v) and not COMBO.search(v)]
        if not ctrl_vals or not drug_vals:
            continue
        ctrl = vals.index[vals.isin(ctrl_vals)].tolist()
        treat = vals.index[vals.isin(drug_vals)].tolist()
        if len(ctrl) >= 2 and len(treat) >= 2:
            cand = (ctrl, treat, field)
            if best is None or len(treat) + len(ctrl) > len(best[0]) + len(best[1]):
                best = cand
    return best


def _labels_for(samples: pd.DataFrame, gsm: str) -> list[str]:
    """Every name GEO gives a sample that a supplementary table might use as a column.

    Includes the title, the ``description`` field (which often carries an explicit
    "Library name: NC2_1"), and the basename of its supplementary file.
    """
    out = [gsm]
    row = samples.loc[gsm]
    for field in ("title", "description", "source_name_ch1"):
        v = row.get(field)
        if isinstance(v, str) and v.strip():
            out.append(v.strip())
            m = re.match(r"(?:library name|sample name|library)\s*[:=]\s*(.+)$", v.strip(), re.I)
            if m:
                out.append(m.group(1).strip())
    for field in [c for c in samples.columns if c.startswith("supplementary_file")]:
        v = row.get(field)
        if isinstance(v, str) and v.strip() and v.strip().upper() != "NONE":
            base = v.strip().split("/")[-1]
            while re.search(r"\.(gz|bz2|zip|txt|tsv|csv|xlsx?|cel|bam|fastq)$", base, re.I):
                base = re.sub(r"\.[^.]+$", "", base)
            out.append(base)
    return out


def _norm_label(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def match_columns(samples: pd.DataFrame, gsms: list[str], columns: list[str]) -> dict[str, str] | None:
    """Map each GSM to exactly one column of a supplementary table, or None.

    Requires a unique match for every requested sample. Anything ambiguous is
    refused rather than guessed, so a mislabelled arm cannot enter a signature.
    """
    norm_cols: dict[str, list[str]] = {}
    for c in columns:
        norm_cols.setdefault(_norm_label(c), []).append(c)
    out: dict[str, str] = {}
    for gsm in gsms:
        hits: list[str] = []
        for lab in _labels_for(samples, gsm):
            key = _norm_label(lab)
            if not key:
                continue
            if key in norm_cols and len(norm_cols[key]) == 1:
                hits.append(norm_cols[key][0])
                continue
            # Either may embed the other: a column "run_NC2_1_count" contains the
            # label "NC2_1", while a title "Non Treated replicate 1 (NT1)" contains
            # the column "NT1". Require >= 3 characters so short tokens cannot match
            # everything, and a unique hit either way.
            if len(key) >= 3:
                part = [c for k, cs in norm_cols.items() for c in cs
                        if key in k or (len(k) >= 3 and k in key)]
                if len(set(part)) == 1:
                    hits.append(part[0])
        uniq = sorted(set(hits))
        if len(uniq) > 1:
            # A table often carries the same sample twice, as a count and an FPKM/TPM
            # column. Those are not ambiguous labels, just two value types: prefer the
            # counts. Anything still ambiguous is refused.
            counts = [c for c in uniq if re.search(r"count", c, re.I)]
            uniq = counts if len(counts) == 1 else uniq
        if len(uniq) != 1:
            return None
        out[gsm] = uniq[0]
    if len(set(out.values())) != len(out):
        return None
    return out


# GEO expression tables often carry whole GO/KEGG annotation blobs in one field, well
# past the csv module's default 128 KiB limit, which otherwise aborts the read.
csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

PROCESSED_SUPP = re.compile(r"(count|fpkm|tpm|cpm|rpkm|matrix|expression|normali[sz]ed|rma)", re.I)
GENE_COL = re.compile(r"^(gene[_ ]?name|gene[_ ]?symbol|symbol|gene[_ ]?id|geneid|id|ensembl[_ ]?id|"
                      r"gene|feature[_ ]?id|test_id)$", re.I)


def supplementary_table(series_dir: Path, samples: pd.DataFrame, need: list[str] | None = None,
                        why: list[str] | None = None, symbols: set[str] | None = None) -> pd.DataFrame | None:
    """Genes x GSM expression from a processed supplementary table, or None.

    Columns are mapped to samples via GEO's own labels; a table that cannot be
    mapped unambiguously is skipped, not guessed at. Only ``need`` (the samples the
    arms actually use) must map, because a series often publishes one file per cell
    line and no single file covers every sample.
    """
    need = list(need) if need else list(samples.index)
    files = [p for p in sorted(series_dir.iterdir())
             if p.is_file() and "series_matrix" not in p.name and PROCESSED_SUPP.search(p.name)
             and p.suffix in (".gz", ".txt", ".tsv", ".csv")]
    if not files:
        if why is not None:
            why.append("no expression table in the series matrix and no processed supplement")
        return None
    for p in files:
        for enc in ("utf-8", "utf-16"):
            try:
                head = pd.read_csv(p, sep=None, engine="python", encoding=enc, nrows=3)
            except Exception:
                continue
            cols = [str(c).strip() for c in head.columns]
            mapping = match_columns(samples, need, cols)
            if mapping is None:
                continue
            try:
                df = pd.read_csv(p, sep=None, engine="python", encoding=enc)
            except Exception:
                continue
            df.columns = [str(c).strip() for c in df.columns]
            # Choose the identifier column by CONTENT, not by name: a table often has
            # both "Gene id" (Ensembl) and "Gene name" (symbol), and only symbols can be
            # compared across series.
            cands = [c for c in df.columns if GENE_COL.match(str(c))] or [df.columns[0]]
            if symbols is not None:
                scored = [(df[c].astype(str).str.upper().str.strip().isin(symbols).mean(), c)
                          for c in cands]
                gene_col = max(scored)[1]
            else:
                gene_col = cands[0]
            expr = df[list(mapping.values())].apply(pd.to_numeric, errors="coerce")
            idx = df[gene_col].astype(str)
            if "_" in idx.iloc[0] and idx.str.startswith(("ENSG", "ENST")).mean() > 0.5:
                idx = idx.str.split("_").str[-1]     # rows like ENSG..._SYMBOL
            expr.index = idx.str.upper().str.strip()
            expr = expr[(expr.index != "") & ~expr.index.isin(["NAN", "NA", "NONE", "-"])]
            expr = expr.groupby(level=0).sum()
            expr.columns = list(mapping.keys())      # back to GSM ids
            if why is not None:
                why.append(f"used supplementary table {p.name}")
            return expr
    if why is not None:
        why.append("processed supplement present but its columns could not be matched to samples")
    return None


def load_symbol_reference(raw_dir: Path) -> set[str] | None:
    """Known human gene symbols, from the NCBI map `fetch-geo --gene-info` writes."""
    p = raw_dir / "ensembl_to_symbol.tsv"
    if not p.exists():
        return None
    s = pd.read_csv(p, sep="\t", header=None, usecols=[1], dtype=str)[1]
    return set(s.str.upper())


def series_signature(gse: str, drug: str, raw_dir: Path, min_genes: int = 3000,
                     symbols: set[str] | None = None, min_symbol_frac: float = 0.3,
                     why: list[str] | None = None):
    """One drug signature from a GEO series.

    Expression comes from the series matrix where it carries a table (microarray),
    otherwise from a processed supplementary table whose columns are matched to
    samples through GEO's own labels (:func:`match_columns`), which refuses anything
    ambiguous rather than guessing.
    """
    fetch_geo_series(gse, raw_dir, suppl=True, log=lambda *a, **k: None)
    mats = sorted((raw_dir / gse).glob("*series_matrix.txt.gz"))
    for m in mats:
        expr, samples = read_series_matrix(m)
        arms = detect_arms(samples, drug)
        if not arms:
            if why is not None:
                why.append("no control/treated arms in the sample annotations")
            continue
        ctrl, treat, field = arms
        if expr.shape[0] < 500:
            expr = supplementary_table(raw_dir / gse, samples, need=ctrl + treat, why=why,
                                       symbols=symbols)
            if expr is None:
                continue
        pm_path = (raw_dir / gse / "probe_map.tsv")
        if not pm_path.exists():
            cands = sorted((raw_dir / gse).glob("probe_map_*.tsv"))
            pm_path = cands[0] if cands else None
        if pm_path is not None and pm_path.exists():
            pm = pd.read_csv(pm_path, sep="\t", header=None, index_col=0, dtype=str)[1]
            expr = collapse_probes(expr, pm)
        else:
            expr.index = expr.index.astype(str).str.upper()
            expr = expr.groupby(level=0).mean()
        expr = maybe_log2(expr)
        expr = expr.dropna(how="all")
        if expr.shape[0] < min_genes:
            if why is not None:
                why.append(f"only {expr.shape[0]} rows after probe collapse")
            continue
        # The rows must be gene symbols. Without a platform probe map they stay as
        # probe ids, which share no identifiers with other series, so the consensus
        # would silently have nothing to average. Check against known symbols.
        if symbols is not None:
            hit = expr.index.astype(str).str.upper().isin(symbols)
            if hit.mean() < min_symbol_frac or int(hit.sum()) < min_genes:
                if why is not None:
                    why.append(f"rows are not gene symbols ({hit.mean():.0%} matched); no probe map")
                continue
            expr = expr[hit]
        meta = {"npi_id": f"{drug}|{gse}", "contrast": f"{drug} vs control ({field}), unpaired",
                "provenance": "geo_drug_treatment", "source_accessions": [gse]}
        sig = build_signature(expr[ctrl], expr[treat], meta, paired=False, min_n=2)
        sig.kind = DRUG
        sig.modality = "compound"
        sig.meta = dict(sig.meta or {}) | {"drug_id": drug, "field": field,
                                           "n_control": len(ctrl), "n_treated": len(treat)}
        return sig
    return None


def consensus(sigs: list, drug: str, min_series: int = 2, min_shared: int = 2000):
    """Per-gene median over series signatures, re-standardised."""
    if len(sigs) < min_series:
        return None
    frame = pd.concat([s.as_series() for s in sigs], axis=1, join="inner")
    if len(frame) < min_shared:
        return None
    med = standardize(frame.median(axis=1))
    out = sigs[0]
    base = type(out)(
        sig_id=drug, kind=DRUG, genes=list(med.index), z=med.to_numpy(),
        provenance="geo_drug_consensus", modality="compound", tissue="multi_cell_line",
        species="human", contrast=f"median of {len(sigs)} independent GEO series",
        quality_flag="ok", source_accessions=sorted({a for s in sigs for a in s.source_accessions}),
        meta={"drug_id": drug, "n_series": len(sigs),
              "series": [s.meta.get("field") and s.source_accessions[0] for s in sigs],
              "n_genes": int(len(med))},
    )
    return base


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--drugs", nargs="*", default=MISSING_FROM_LINCS,
                    help="default: the drugs LINCS Phase II does not carry")
    ap.add_argument("--max-per-drug", type=int, default=60)
    ap.add_argument("--max-series-used", type=int, default=6)
    ap.add_argument("--min-score", type=int, default=5)
    ap.add_argument("--array-only", action="store_true", default=True,
                    help="microarray series only: their series matrix carries the expression table, so "
                         "sample columns are GSMs that GEO's annotation labels. RNA-seq series put "
                         "counts in supplementary files under arbitrary column names.")
    ap.add_argument("--allow-rnaseq", dest="array_only", action="store_false")
    ap.add_argument("--raw-dir", default="data/raw")
    ap.add_argument("--validate-against", default="data/processed/signatures/lincs_all.parquet")
    ap.add_argument("--out", default="data/processed/signatures/drugs_consensus.parquet")
    ap.add_argument("--report", default="out/drug_consensus")
    a = ap.parse_args(argv)
    raw = Path(a.raw_dir)
    rep = Path(a.report)
    rep.mkdir(parents=True, exist_ok=True)

    symbols = load_symbol_reference(raw)
    if symbols is None:
        print("warning: no data/raw/ensembl_to_symbol.tsv; run npi-pharma fetch-geo --gene-info", file=sys.stderr)
    per_series, consensuses, log = {}, [], []
    for drug in a.drugs:
        types = ('"expression profiling by array"[gdstyp]' if a.array_only else
                 '("expression profiling by array"[gdstyp] OR '
                 '"expression profiling by high throughput sequencing"[gdstyp])')
        term = f'"{drug}" AND "Homo sapiens"[orgn] AND gse[etyp] AND {types}'
        uids = esearch(term, a.max_per_drug)
        recs = []
        for i in range(0, len(uids), 100):
            recs += esummary(uids[i:i + 100])
            time.sleep(0.5)
        recs = [r for r in recs if score(r, drug)[0] >= a.min_score]
        recs.sort(key=lambda r: -score(r, drug)[0])
        print(f"{drug}: {len(recs)} candidate series", file=sys.stderr)
        got = []
        for r in recs:
            if len(got) >= a.max_series_used:
                break
            acc = r.get("accession")
            why: list[str] = []
            try:
                sig = series_signature(acc, drug, raw, symbols=symbols, why=why)
            except Exception as e:
                log.append({"drug": drug, "accession": acc, "outcome": f"error: {type(e).__name__}: {e}"})
                continue
            if sig is None:
                log.append({"drug": drug, "accession": acc,
                            "outcome": why[0] if why else "no matrix file"})
                continue
            got.append(sig)
            log.append({"drug": drug, "accession": acc, "outcome": "used",
                        "field": sig.meta.get("field"), "n_control": sig.meta.get("n_control"),
                        "n_treated": sig.meta.get("n_treated"), "n_genes": len(sig.genes)})
            print(f"  {acc}: {sig.meta['n_control']}v{sig.meta['n_treated']} "
                  f"on '{sig.meta['field']}', {len(sig.genes)} genes", file=sys.stderr)
        per_series[drug] = got
        c = consensus(got, drug)
        if c is not None:
            consensuses.append(c)

    pd.DataFrame(log).to_csv(rep / "series_log.tsv", sep="\t", index=False)
    if per_series:
        flat = [s for v in per_series.values() for s in v]
        if flat:
            save_signatures(flat, str(Path(a.out).with_name("drugs_geo_per_series.parquet")))
    if consensuses:
        save_signatures(consensuses, a.out)

    # QC: cross-series agreement within a drug, and agreement with LINCS where both exist
    qc = {}
    for drug, sigs in per_series.items():
        if len(sigs) < 2:
            continue
        f = pd.concat([s.as_series() for s in sigs], axis=1, join="inner")
        z = f.apply(lambda c: (c - c.mean()) / c.std(ddof=0))
        C = (z.T @ z) / len(z)
        iu = np.triu_indices_from(C, 1)
        qc[drug] = {"n_series": len(sigs), "shared_genes": int(len(z)),
                    "median_cross_series_cosine": float(np.median(C.to_numpy()[iu]))}
    lincs_path = Path(a.validate_against)
    if lincs_path.exists() and consensuses:
        lincs = {s.sig_id.lower(): s.as_series() for s in load_signatures(lincs_path)}
        for c in consensuses:
            key = c.sig_id.lower()
            key = {"mitomycin c": "mitomycin-c"}.get(key, key)
            if key in lincs:
                x, y = c.as_series(), lincs[key]
                g = x.index.intersection(y.index)
                if len(g) >= 300:
                    xs = (x[g] - x[g].mean()) / x[g].std(ddof=0)
                    ys = (y[g] - y[g].mean()) / y[g].std(ddof=0)
                    qc.setdefault(c.sig_id, {})["vs_lincs_consensus_cosine"] = float((xs @ ys) / len(g))
                    qc[c.sig_id]["vs_lincs_shared_genes"] = int(len(g))
    (rep / "qc.json").write_text(json.dumps(qc, indent=2, default=float))

    print(f"\nconsensus signatures built: {len(consensuses)} "
          f"({', '.join(c.sig_id for c in consensuses)})")
    for drug, v in qc.items():
        bits = ", ".join(f"{k}={v[k]:.3f}" if isinstance(v[k], float) else f"{k}={v[k]}" for k in v)
        print(f"  {drug}: {bits}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
