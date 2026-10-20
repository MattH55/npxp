"""Pediatric context for mechanism-hypothesis genes (lookups, not evidence of interaction).

For each gene and each pediatric cancer group with enough samples:
  - putative LoF SNV/indel fraction (WGS/WXS tumours; panels excluded so a
    gene outside a panel cannot read as "never mutated")
  - deep-deletion fraction (consensus CNV, same DNA samples)
  - fusion fraction (putative-oncogenic fusions, RNA tumours)
  - expression median / IQR (log2 TPM+1) and the group median's percentile
    across groups
plus the gene's FDA PMTL designation from OpenPedCan's mapping.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .openpedcan import LOF_CLASSES

PEDIATRIC_COHORTS = {"PBTA", "TARGET", "GMKF", "DGD", "PPTC", "Maris"}
# Only for this category is the lesion gene a loss-of-function state.
LOF_CATEGORY = "induced_hr_deficiency"
NOTE = "Pediatric context: annotation from OpenPedCan, not evidence of an NPI x drug interaction."


def dna_samples(hist: pd.DataFrame, independent_dna: pd.DataFrame) -> pd.DataFrame:
    """One WGS/WXS primary tumour per patient (pediatric cohorts, solid/blood/marrow, no models)."""
    keep = set(independent_dna.loc[independent_dna["experimental_strategy"].isin(["WGS", "WXS"]),
                                   "Kids_First_Biospecimen_ID"])
    h = hist[hist["Kids_First_Biospecimen_ID"].isin(keep) & hist["cohort"].isin(PEDIATRIC_COHORTS)]
    h = h[~h["composition"].isin({"Derived Cell Line", "Patient Derived Xenograft"}) & h["cancer_group"].notna()]
    return h[["Kids_First_Biospecimen_ID", "cancer_group", "cohort", "experimental_strategy"]]


def read_maf_genes(maf_path: str | Path, genes: set[str], chunksize: int = 500_000) -> pd.DataFrame:
    cols = ["Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode", "HGVSp_Short"]
    parts = []
    for chunk in pd.read_csv(maf_path, sep="\t", usecols=cols, dtype=str, chunksize=chunksize, comment="#"):
        parts.append(chunk[chunk["Hugo_Symbol"].isin(genes)])
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=cols)


def read_cnv_genes(cnv_path: str | Path, genes: set[str], chunksize: int = 2_000_000) -> tuple[pd.DataFrame, set[str]]:
    """Deep-deletion calls for ``genes`` and the set of biospecimens present in the CNV file."""
    cols = ["biospecimen_id", "status", "gene_symbol"]
    parts, seen = [], set()
    for chunk in pd.read_csv(cnv_path, sep="\t", usecols=cols, dtype=str, chunksize=chunksize):
        seen.update(chunk["biospecimen_id"].unique())
        parts.append(chunk[chunk["gene_symbol"].isin(genes) & (chunk["status"] == "deep deletion")])
    return pd.concat(parts, ignore_index=True), seen


def fusion_genes(fusions: pd.DataFrame) -> pd.DataFrame:
    """Long table (Sample, gene) of genes partnered in a putative-oncogenic fusion."""
    gcols = [c for c in ("Gene1A", "Gene1B", "Gene2A", "Gene2B") if c in fusions]
    long = fusions.melt(id_vars=["Sample"], value_vars=gcols, value_name="gene").dropna(subset=["gene"])
    return long[["Sample", "gene"]].drop_duplicates()


def _frac(hits: pd.DataFrame, sample_col: str, gene_col: str, groups: dict[str, set[str]], genes: list[str]) -> dict:
    out = {}
    for grp, ss in groups.items():
        sub = hits[hits[sample_col].isin(ss)]
        n = sub.groupby(gene_col)[sample_col].nunique()
        out[grp] = {g: (int(n.get(g, 0)), len(ss)) for g in genes}
    return out


def gene_context(
    genes: list[str],
    hist: pd.DataFrame,
    rna: pd.DataFrame,
    expr: pd.DataFrame,
    independent_dna: pd.DataFrame,
    maf: pd.DataFrame,
    cnv_deep: pd.DataFrame,
    cnv_samples: set[str],
    fusions: pd.DataFrame,
    pmtl: pd.DataFrame,
    min_n: int = 20,
) -> pd.DataFrame:
    """Gene x cancer-group long table. ``rna`` = tumour rows from cohorts.rna_samples."""
    dna = dna_samples(hist, independent_dna)
    mut_groups = {g: set(d["Kids_First_Biospecimen_ID"]) for g, d in dna.groupby("cancer_group") if len(d) >= min_n}
    cnv_groups = {g: s & cnv_samples for g, s in mut_groups.items()}
    cnv_groups = {g: s for g, s in cnv_groups.items() if len(s) >= min_n}
    tum = rna[rna["kind"] == "tumour"]
    rna_groups = {g: set(d["Kids_First_Biospecimen_ID"]) for g, d in tum.groupby("cancer_group") if len(d) >= min_n}

    lof = maf[maf["Variant_Classification"].isin(LOF_CLASSES)]
    f_lof = _frac(lof, "Tumor_Sample_Barcode", "Hugo_Symbol", mut_groups, genes)
    f_del = _frac(cnv_deep, "biospecimen_id", "gene_symbol", cnv_groups, genes)
    f_fus = _frac(fusion_genes(fusions), "Sample", "gene", rna_groups, genes)

    present = [g for g in genes if g in expr.index]
    med = {grp: expr.loc[present, sorted(ss)].median(axis=1) for grp, ss in rna_groups.items()}
    q1 = {grp: expr.loc[present, sorted(ss)].quantile(0.25, axis=1) for grp, ss in rna_groups.items()}
    q3 = {grp: expr.loc[present, sorted(ss)].quantile(0.75, axis=1) for grp, ss in rna_groups.items()}
    med_df = pd.DataFrame(med)
    pct = med_df.rank(axis=1, pct=True) * 100 if len(med_df.columns) else med_df

    pm = pmtl.drop_duplicates("gene_symbol").set_index("gene_symbol")["pmtl"]
    rows = []
    for grp in sorted(set(mut_groups) | set(rna_groups)):
        for g in genes:
            r = {"cancer_group": grp, "gene": g, "pmtl": None if pd.isna(pm.get(g)) else pm.get(g)}
            if grp in f_lof:
                k, n = f_lof[grp][g]
                r.update(n_dna=n, n_lof=k, lof_fraction=k / n)
            if grp in f_del:
                k, n = f_del[grp][g]
                r.update(n_cnv=n, n_deep_del=k, deep_del_fraction=k / n)
            if grp in f_fus:
                k, n = f_fus[grp][g]
                r.update(n_rna=n, n_fusion=k, fusion_fraction=k / n)
            if grp in med and g in med_df.index:
                r.update(expr_median=float(med[grp][g]), expr_q1=float(q1[grp][g]), expr_q3=float(q3[grp][g]),
                         expr_group_percentile=float(pct.at[g, grp]))
            rows.append(r)
    return pd.DataFrame(rows)


def row_summary(row: dict, ctx: pd.DataFrame, lof_threshold: float = 0.05, top_k: int = 3) -> dict:
    """Compact per-hypothesis block: where the lesion is already present at baseline,
    and where the drug target is most expressed. Existing row fields are untouched."""
    out = {"note": NOTE}
    lesion = row.get("lesion_gene")
    if lesion and row.get("category", LOF_CATEGORY) != LOF_CATEGORY:
        out["lesion_note"] = ("baseline LoF/deep-deletion frequency is not reported: for this category the "
                              "lesion gene is a pathway or immune target, not a loss-of-function lesion")
    elif lesion:
        c = ctx[ctx["gene"] == lesion].copy()
        c["altered"] = c[["lof_fraction", "deep_del_fraction"]].fillna(0).max(axis=1)
        hit = c[c["altered"] >= lof_threshold].sort_values("altered", ascending=False)
        out["lesion_already_altered"] = [
            {"cancer_group": r.cancer_group, "lof_fraction": _r(r.lof_fraction), "deep_del_fraction": _r(r.deep_del_fraction)}
            for r in hit.itertuples()]
        top = c.dropna(subset=["n_dna"]).sort_values("altered", ascending=False).head(1)
        out["lesion_max_altered"] = None if top.empty else {
            "gene": lesion, "cancer_group": top.iat[0, top.columns.get_loc("cancer_group")],
            "altered_fraction": _r(top["altered"].iat[0]), "n_dna": int(top["n_dna"].iat[0])}
    tgt = row.get("drug_target_gene")
    if tgt:
        c = ctx[(ctx["gene"] == tgt) & ctx["expr_median"].notna()].sort_values("expr_median", ascending=False)
        out["target_expression_top"] = [
            {"cancer_group": r.cancer_group, "median_log2tpm": _r(r.expr_median)} for r in c.head(top_k).itertuples()]
        p = ctx.loc[ctx["gene"] == tgt, "pmtl"].dropna()
        out["target_pmtl"] = p.iat[0] if len(p) else None
    return out


def _r(x, nd: int = 3):
    return None if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), nd)
