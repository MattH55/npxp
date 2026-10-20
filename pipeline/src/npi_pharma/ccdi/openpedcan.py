"""OpenPedCan (Open Pediatric Cancer) release files from the public D3b bucket.

OpenPedCan harmonises TARGET, PBTA/CBTN, GMKF (Kids First) and other pediatric
cohorts with one pipeline (RSEM TPM, consensus SNV/CNV, annoFuse fusions) and
ships the FDA Pediatric Molecular Target List mapping. Open access, CC-BY-4.0.
The release's md5sum.txt is the source of truth for file names and checksums.
"""

from __future__ import annotations

import hashlib
import shutil
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from ..genes import normalize_symbol
from . import manifest

BUCKET = "https://s3.amazonaws.com/d3b-openaccess-us-east-1-prd-pbta/open-targets"
DEFAULT_RELEASE = "v15"  # last complete release; v16 is an incremental drop without histologies
CITATION = ("OpenPedCan: Geng Z. et al., 'The Open Pediatric Cancer Project', bioRxiv 2024; "
            "data CC-BY-4.0, https://github.com/d3b-center/OpenPedCan-analysis")

CORE_FILES = [
    "histologies.tsv",
    "ensg-hugo-pmtl-mapping.tsv",
    "independent-specimens.rnaseqpanel.primary.tsv",
    "independent-specimens.wgswxspanel.primary.prefer.wgs.tsv",
    "gene-expression-rsem-tpm-collapsed.rds",
    "fusion-putative-oncogenic.tsv",
    "consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only.tsv.gz",
    "release-notes.md",
]
MAF = "snv-consensus-plus-hotspots.maf.tsv.gz"

# Variant classes counted as putative loss of function.
LOF_CLASSES = {"Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Nonstop_Mutation",
               "Translation_Start_Site"}


def read_md5(text: str) -> dict[str, str]:
    out = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 2:
            out[parts[-1]] = parts[0]
    return out


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fetch(raw_dir: str | Path, files: list[str] | None = None, release: str = DEFAULT_RELEASE,
          overwrite: bool = False, log=print) -> dict[str, Path]:
    """Download release files, verify md5 against the release's md5sum.txt, write MANIFEST.json."""
    out = Path(raw_dir) / "openpedcan" / release
    out.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(f"{BUCKET}/{release}/md5sum.txt", timeout=60) as r:
        md5_text = r.read().decode()
    (out / "md5sum.txt").write_text(md5_text)
    md5s = read_md5(md5_text)
    got, entries = {}, {}
    for name in files or CORE_FILES:
        dest = out / name
        url = f"{BUCKET}/{release}/{name}"
        expect = md5s.get(name)
        if dest.exists() and not overwrite and (expect is None or _md5(dest) == expect):
            log(f"  have {name}")
        else:
            tmp = dest.with_suffix(dest.suffix + ".part")
            with urllib.request.urlopen(url, timeout=300) as r, open(tmp, "wb") as fh:
                shutil.copyfileobj(r, fh, length=1 << 22)
            if expect is not None and _md5(tmp) != expect:
                tmp.unlink()
                raise ValueError(f"md5 mismatch for {name}")
            tmp.replace(dest)
            log(f"  got {name} ({dest.stat().st_size / 1e6:.1f} MB)")
        got[name] = dest
        entries[name] = {"url": url, "md5": expect, "sha256": manifest.file_digest(dest),
                         "size": dest.stat().st_size, "retrieved": manifest.now_utc()}
    manifest.record(out, {"source": "OpenPedCan", "release": release, "access": "open",
                          "license": "CC-BY-4.0", "citation": CITATION, "bucket": BUCKET}, entries)
    return got


# ---------------------------------------------------------------- readers

def read_histologies(path: str | Path) -> pd.DataFrame:
    h = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    h["age_at_diagnosis_years"] = pd.to_numeric(h.get("age_at_diagnosis_days"), errors="coerce") / 365.25
    return h


def read_tpm_rds(path: str | Path, samples: list[str] | None = None) -> pd.DataFrame:
    """Gene-symbol x sample TPM from the collapsed RSEM .rds (an R data.frame)."""
    import pyreadr

    df = next(iter(pyreadr.read_r(str(path)).values()))
    if "gene_id" in df.columns:  # some releases carry the symbol as a column
        df = df.set_index("gene_id")
    if samples is not None:
        df = df[[s for s in samples if s in df.columns]]
    df.index = [normalize_symbol(g) for g in df.index]
    return df.astype(np.float32)


def primary_rna_samples(hist: pd.DataFrame, independent: pd.DataFrame | None = None) -> pd.DataFrame:
    """Tumour RNA-seq biospecimens, one per patient when the independent-specimens list is given."""
    h = hist[(hist["sample_type"] == "Tumor") & (hist["experimental_strategy"] == "RNA-Seq")]
    if independent is not None:
        h = h[h["Kids_First_Biospecimen_ID"].isin(set(independent["Kids_First_Biospecimen_ID"]))]
    return h


def cohort_table(hist: pd.DataFrame, min_n: int = 20, key: str = "cancer_group") -> pd.DataFrame:
    """Counts per cohort / cancer group, flagging groups below ``min_n``."""
    t = hist.groupby([key, "cohort"]).size().rename("n").reset_index()
    tot = hist.groupby(key).size().rename("n_total")
    t = t.merge(tot, left_on=key, right_index=True)
    t["enough"] = t["n_total"] >= min_n
    return t.sort_values("n_total", ascending=False)


def lof_frequency(maf: pd.DataFrame, genes: list[str], samples_by_group: dict[str, set[str]]) -> pd.DataFrame:
    """Fraction of DNA tumour samples per group with a putative LoF SNV/indel in each gene.

    ``maf`` needs Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode.
    """
    m = maf[maf["Variant_Classification"].isin(LOF_CLASSES) & maf["Hugo_Symbol"].isin(genes)]
    rows = []
    for grp, samples in samples_by_group.items():
        sub = m[m["Tumor_Sample_Barcode"].isin(samples)]
        hit = sub.groupby("Hugo_Symbol")["Tumor_Sample_Barcode"].nunique()
        for g in genes:
            k = int(hit.get(g, 0))
            rows.append({"group": grp, "gene": g, "n_samples": len(samples), "n_lof": k,
                         "lof_fraction": k / len(samples) if samples else np.nan})
    return pd.DataFrame(rows)
