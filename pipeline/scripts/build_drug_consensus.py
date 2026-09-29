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


# "_" is a word character, so \b does not separate "DMSO_R1"; match on non-alphanumeric instead.
_REPLICATE = re.compile(r"(?<![A-Za-z0-9])(?:r|rep|replicate)[ _-]?(\d+)(?![A-Za-z0-9])", re.I)


def _tokens(s: str) -> list[str]:
    """Label to comparable tokens, splitting run-together words and digit runs.

    ``CAMA-1_DMSO_1`` and ``SNB19ACT`` both become word/number tokens, so a column
    written without separators can still be compared with a spaced GEO title.
    ``R1``/``Rep1`` are rewritten to ``replicate 1`` first, because a table writes
    the replicate number as ``_R1`` where the title spells it out.
    """
    s = _REPLICATE.sub(r"replicate \1", str(s))
    s = re.sub(r"(?<=[a-zA-Z])(?=[0-9])|(?<=[0-9])(?=[a-zA-Z])", " ", s)
    s = re.sub(r"(?<=[a-z])(?=[A-Z])", " ", s)
    return [t for t in re.split(r"[^A-Za-z0-9]+", s.lower()) if t]


def _token_score(column: str, label: str) -> float:
    """Fraction of a column's tokens accounted for by a sample label's tokens.

    A prefix counts partially: tables abbreviate (``Veh`` for ``Vehicle``,
    ``Sta`` for ``Stattic``) but do not rename. Prefixes must be >= 3 characters
    so that a stray letter cannot match a word.
    """
    ct, lt = _tokens(column), set(_tokens(label))
    if not ct:
        return 0.0
    total = 0.0
    for t in ct:
        if t in lt:
            total += 1.0
        elif len(t) >= 3 and any(u.startswith(t) or (len(u) >= 3 and t.startswith(u)) for u in lt):
            total += 0.6
    return total / len(ct)


def match_columns_by_assignment(samples: pd.DataFrame, gsms: list[str], columns: list[str],
                                min_score: float = 0.8) -> dict[str, str] | None:
    """Map GSMs to columns as a one-to-one assignment, or None if it is not decisive.

    Many RNA-seq series publish a count table whose columns are *arm labels*
    (``CAMA-1_DMSO_1``, ``Veh_R``) rather than GSM names, so no single label
    matches a column outright. Solving all the samples at once is what makes that
    tractable: a column has to beat every other sample as well as every other
    column, and the control arm is often identified only by which column is left
    over. Three conditions must all hold, or the table is refused:

      * every sample scores at least ``min_score`` against its assigned column;
      * the assignment is one-to-one;
      * forbidding any assigned pair makes the total strictly worse, so no second
        arrangement scores as well. This is the guard against guessing -- a table
        whose labels are genuinely ambiguous has a tie, and a tie is a refusal.
    """
    from scipy.optimize import linear_sum_assignment

    # Its whole justification is solving the samples jointly -- a column must beat
    # every other sample as well as every other column. With one sample there is no
    # such constraint, so a single-sample request is refused rather than guessed.
    if len(gsms) < 2 or len(columns) < len(gsms):
        return None
    # A table often publishes each sample twice, as counts and as TPM/FPKM. Those are
    # two value types of one sample, not two candidate labels, and they tie exactly on
    # token score -- which would refuse the whole table. Keep the counts.
    # The suffix is also not part of the sample's name, so scoring uses the stem.
    stems: dict[str, list[str]] = {}
    for c in columns:
        stems.setdefault(re.sub(r"[ _.-]?(counts?|tpms?|fpkms?|cpms?|rpkms?)$", "", c,
                                flags=re.I), []).append(c)
    keys = list(stems)
    pick = [cs[0] if len(cs) == 1 else
            next((c for c in cs if re.search(r"counts?$", c, re.I)), cs[0])
            for cs in stems.values()]
    if len(keys) < len(gsms):
        return None
    labels = [" ".join(_labels_for(samples, g)) for g in gsms]
    S = np.array([[_token_score(k, lab) for k in keys] for lab in labels])

    def solve(mask: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
        cost = np.where(mask, -S, 1e6)
        r, c = linear_sum_assignment(cost)
        return float(S[r, c][mask[r, c]].sum() - 1e6 * (~mask[r, c]).sum()), r, c

    mask = np.ones_like(S, dtype=bool)
    best, rows, cols = solve(mask)
    if S[rows, cols].min() < min_score:
        return None
    for r, c in zip(rows, cols):
        m = mask.copy()
        m[r, c] = False
        alt, _, _ = solve(m)
        if alt >= best - 1e-9:
            return None  # a second arrangement does as well: ambiguous, so refuse
    return {gsms[r]: pick[c] for r, c in zip(rows, cols)}


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
            return match_columns_by_assignment(samples, gsms, columns)
        out[gsm] = uniq[0]
    if len(set(out.values())) != len(out):
        return match_columns_by_assignment(samples, gsms, columns)
    return out


# GEO expression tables often carry whole GO/KEGG annotation blobs in one field, well
# past the csv module's default 128 KiB limit, which otherwise aborts the read.
csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

def to_numeric_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Coerce to numbers, retrying a column that is written with comma decimals.

    GEO carries tables exported under a European locale: GSE304295's TPM file is
    semicolon-separated with values like ``2,74724752907299``. Plain coercion turns
    every one of those into NaN, and because the caller then sums by gene -- and
    ``groupby().sum()`` of all-NaN is 0.0, not NaN -- the failure arrives as a table
    of zeros instead of an error.
    """
    out = {}
    for c in df.columns:
        raw = df[c]
        num = pd.to_numeric(raw, errors="coerce")
        nonblank = raw.astype(str).str.strip().replace({"": None}).notna()
        lost = nonblank & num.isna()
        if nonblank.sum() and lost.sum() / nonblank.sum() > 0.5:
            retry = pd.to_numeric(
                raw.astype(str).str.strip().str.replace(".", "", regex=False)
                   .str.replace(",", ".", regex=False), errors="coerce")
            if retry.notna().sum() > num.notna().sum():
                num = retry
        out[c] = num
    return pd.DataFrame(out, index=df.index)


def is_degenerate(expr: pd.DataFrame, min_varying_frac: float = 0.01) -> bool:
    """Is this expression table constant, so it can carry no signal?

    A table that failed to parse arrives as all zeros rather than as an error (see
    :func:`to_numeric_frame`). Standardising it later yields 0/0, so one such series
    turns its drug's whole cross-series median into NaN -- which then reads as a
    measured-low agreement rather than as a broken input.
    """
    if expr.empty or expr.shape[1] < 2:
        return True
    a = expr.to_numpy(float)
    if not np.isfinite(a).any():
        return True
    varying = np.nanstd(a, axis=1) > 0
    return bool(varying.mean() < min_varying_frac)


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
            # The identifier may not be a column at all. A table whose first column has
            # no header -- GSE304295's raw counts -- is read with the gene names as the
            # INDEX, so scoring only columns picks a sample and the symbol gate then
            # rejects the series at 0% matched. Score the index on the same footing.
            index_ids = pd.Series(df.index.astype(str), index=df.index)
            if symbols is not None:
                # On a tie an explicitly named column beats the index, hence the 1/0.
                scored = [(df[c].astype(str).str.upper().str.strip().isin(symbols).mean(), 1, c)
                          for c in cands]
                if not isinstance(df.index, pd.RangeIndex):
                    scored.append((index_ids.str.upper().str.strip().isin(symbols).mean(),
                                   0, None))
                gene_col = max(scored)[2]
            else:
                gene_col = cands[0]
            expr = to_numeric_frame(df[list(mapping.values())])
            idx = index_ids.astype(str) if gene_col is None else df[gene_col].astype(str)
            if "_" in idx.iloc[0] and idx.str.startswith(("ENSG", "ENST")).mean() > 0.5:
                idx = idx.str.split("_").str[-1]     # rows like ENSG..._SYMBOL
            elif idx.str.startswith("ENSG").mean() > 0.5:
                # Bare Ensembl IDs, often carrying a version ("ENSG00000000003.15").
                # Left untranslated the series shares no gene with any other, so it
                # silently drops out of every comparison rather than failing loudly.
                emap = load_ensembl_map(series_dir.parent)
                if emap is not None:
                    idx = idx.str.split(".").str[0].map(emap).fillna(idx)
            expr.index = idx.str.upper().str.strip()
            expr = expr[(expr.index != "") & ~expr.index.isin(["NAN", "NA", "NONE", "-"])]
            expr = expr.groupby(level=0).sum()
            expr.columns = list(mapping.keys())      # back to GSM ids
            # groupby().sum() turns an all-NaN group into 0.0, so a table that failed to
            # parse comes back as a silent field of zeros rather than as an error. That
            # produced a completely flat oxaliplatin signature (GSE304295) which passed
            # every downstream gate and turned its drug's cross-series median into NaN.
            if is_degenerate(expr):
                if why is not None:
                    why.append(f"{p.name} parsed to a constant table (no usable numbers)")
                continue
            if why is not None:
                why.append(f"used supplementary table {p.name}")
            return expr
    if why is not None:
        why.append("processed supplement present but its columns could not be matched to samples")
    return None


def load_ensembl_map(raw_dir: Path) -> "pd.Series | None":
    """Ensembl gene ID -> symbol, from the same NCBI map, keyed without a version."""
    p = raw_dir / "ensembl_to_symbol.tsv"
    if not p.exists():
        return None
    m = pd.read_csv(p, sep="\t", header=None, names=["ensembl", "symbol"], dtype=str)
    return m.drop_duplicates("ensembl").set_index("ensembl")["symbol"]


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
        # A flat series standardises to 0/0, so every pair involving it is NaN and a
        # plain median of the pairs is NaN too -- indistinguishable downstream from a
        # genuinely low agreement. Drop those pairs and say how many were computable.
        z = f.apply(lambda c: (c - c.mean()) / c.std(ddof=0))
        C = (z.T @ z) / len(z)
        iu = np.triu_indices_from(C, 1)
        pairs = C.to_numpy()[iu]
        ok = pairs[np.isfinite(pairs)]
        qc[drug] = {"n_series": len(sigs), "shared_genes": int(len(z)),
                    "median_cross_series_cosine": float(np.median(ok)) if len(ok) else None,
                    "n_pairs": int(len(pairs)), "n_pairs_computable": int(len(ok))}
        if len(ok) < len(pairs):
            qc[drug]["flat_series"] = [s.sig_id for s in sigs
                                       if not np.isfinite(np.nanstd(s.z)) or np.nanstd(s.z) == 0]
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
