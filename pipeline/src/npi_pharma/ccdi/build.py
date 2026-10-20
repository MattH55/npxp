"""build-ccdi steps: cohorts, concordance, annotate (offline, from fetched OpenPedCan files)."""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path

import pandas as pd

from . import annotate, cohorts, manifest
from .openpedcan import MAF, read_histologies


def _slug(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", s.lower()).strip("_")


class Inputs:
    def __init__(self, raw_dir: str | Path, processed: str | Path, release: str = "v15"):
        self.raw = Path(raw_dir) / "openpedcan" / release
        self.processed = Path(processed)
        self.release = release
        self.manifest = manifest.load(self.raw)
        if not self.manifest:
            raise FileNotFoundError(f"{self.raw}/MANIFEST.json missing; run fetch-openpedcan first")
        self._hist = self._expr = None

    @property
    def hist(self) -> pd.DataFrame:
        if self._hist is None:
            self._hist = read_histologies(self.raw / "histologies.tsv")
        return self._hist

    @property
    def expr(self) -> pd.DataFrame:
        if self._expr is None:
            pq = self.processed / f"openpedcan_{self.release}_log2tpm.parquet"
            if not pq.exists():
                cohorts.tpm_to_parquet(self.raw / "gene-expression-rsem-tpm-collapsed.rds", pq)
            self._expr = pd.read_parquet(pq)
        return self._expr

    def rna(self) -> pd.DataFrame:
        ind = pd.read_csv(self.raw / "independent-specimens.rnaseqpanel.primary.tsv", sep="\t")
        return cohorts.rna_samples(self.hist, self.expr.columns, ind)

    def provenance(self) -> dict:
        return {"source": "OpenPedCan", "release": self.release, "license": self.manifest.get("license"),
                "citation": self.manifest.get("citation"), "retrieved": self.manifest.get("updated")}


def step_cohorts(inp: Inputs, dmap_path: str | Path, min_tumours: int = 20, log=print) -> dict:
    rna = inp.rna()
    counts = cohorts.cohort_counts(rna)
    dmap = cohorts.load_diagnosis_map(dmap_path)
    mapped, unmapped = cohorts.map_groups(counts.index.tolist(), dmap)
    out = inp.processed / "cohorts"
    out.mkdir(parents=True, exist_ok=True)
    # Cohorts are stratified by RNA library (poly-A vs stranded): a z-score
    # against a mixed cohort mostly measures library chemistry.
    written, strata = {}, []
    tum = rna[rna["kind"] == "tumour"].assign(RNA_library=lambda d: d["RNA_library"].fillna("unknown"))
    for (g, lib), d in tum.groupby(["cancer_group", "RNA_library"]):
        rec = {"cancer_group": g, "RNA_library": lib, "n_tumours": len(d),
               "cohorts": ",".join(sorted(d["cohort"].unique())), "oncotree_code": mapped.get(g, ""), "cohort_file": ""}
        if len(d) >= min_tumours:
            p = out / f"{_slug(g)}__{_slug(lib)}.parquet"
            inp.expr[sorted(d["Kids_First_Biospecimen_ID"])].to_parquet(p)
            written[f"{g} [{lib}]"] = rec["cohort_file"] = str(p)
        strata.append(rec)
    st = pd.DataFrame(strata).sort_values("n_tumours", ascending=False)
    st.to_csv(inp.processed / "cohort_strata.tsv", sep="\t", index=False)
    counts = counts.assign(oncotree_code=[mapped.get(g, "") for g in counts.index])
    counts.to_csv(inp.processed / "cohort_counts.tsv", sep="\t")
    rna.to_csv(inp.processed / "rna_samples.tsv", sep="\t", index=False)
    log(f"{len(written)} tumour cohorts (group x RNA library, >= {min_tumours}) -> {out}")
    log(f"unmapped to OncoTree ({len(unmapped)}; reported, not dropped): {', '.join(unmapped[:15])}"
        + (" ..." if len(unmapped) > 15 else ""))
    return {"cohorts": written, "unmapped": unmapped}


def step_concordance(inp: Inputs, min_tumours: int = 20, n_var: int = 1000, log=print) -> pd.DataFrame:
    rna = inp.rna()
    rna = rna[(rna["kind"] == "tumour") | rna["Kids_First_Biospecimen_ID"].isin(inp.expr.columns)]
    res = cohorts.concordance(inp.expr, rna, min_tumours=min_tumours, n_var=n_var)
    p = inp.processed / "model_concordance.tsv"
    res.to_csv(p, sep="\t", index=False, float_format="%.4f")
    log(f"{len(res)} models -> {p}")
    log(res["flag"].value_counts().to_string())
    return res


def mechanism_genes(rows: list[dict]) -> list[str]:
    g = set()
    for r in rows:
        for k in ("lesion_gene", "sl_partner_gene", "drug_target_gene"):
            if r.get(k):
                g.add(str(r[k]).upper())
    return sorted(g)


def step_annotate(inp: Inputs, hypotheses: str | Path, out_json: str | Path, min_n: int = 20, log=print) -> dict:
    doc = json.loads(Path(hypotheses).read_text())
    # Re-annotating an already-annotated file replaces the old block (idempotent).
    doc.pop("pediatric_context", None)
    rows = [{k: v for k, v in r.items() if k != "pediatric_context"} for r in doc["rows"]]
    doc["rows"] = rows
    genes = mechanism_genes(rows)
    log(f"{len(genes)} genes across {len(rows)} hypothesis rows")
    # Streaming the 2 GB MAF / 0.5 GB CNV takes minutes; cache the gene subset.
    key = hashlib.sha1(",".join(genes).encode()).hexdigest()[:10]
    cache = inp.processed / "cache"
    cache.mkdir(parents=True, exist_ok=True)
    maf_c, cnv_c, seen_c = cache / f"maf_{key}.parquet", cache / f"cnv_{key}.parquet", cache / f"cnv_samples_{key}.txt"
    if maf_c.exists() and cnv_c.exists() and seen_c.exists():
        maf, cnv = pd.read_parquet(maf_c), pd.read_parquet(cnv_c)
        cnv_samples = set(seen_c.read_text().split())
    else:
        maf = annotate.read_maf_genes(inp.raw / MAF, set(genes))
        cnv, cnv_samples = annotate.read_cnv_genes(
            inp.raw / "consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only.tsv.gz", set(genes))
        maf.to_parquet(maf_c)
        cnv.to_parquet(cnv_c)
        seen_c.write_text("\n".join(sorted(cnv_samples)))
    fus = pd.read_csv(inp.raw / "fusion-putative-oncogenic.tsv", sep="\t", dtype=str)
    pmtl = pd.read_csv(inp.raw / "ensg-hugo-pmtl-mapping.tsv", sep="\t", dtype=str)
    ind_dna = pd.read_csv(inp.raw / "independent-specimens.wgswxspanel.primary.prefer.wgs.tsv", sep="\t", dtype=str)
    ctx = annotate.gene_context(genes, inp.hist, inp.rna(), inp.expr, ind_dna, maf, cnv, cnv_samples, fus, pmtl,
                                min_n=min_n)
    ctx.to_csv(inp.processed / "pediatric_gene_context.tsv", sep="\t", index=False, float_format="%.4f")
    before = json.dumps(rows, sort_keys=True)
    for r in rows:
        r["pediatric_context"] = annotate.row_summary(r, ctx)
    # existing fields must be byte-identical
    stripped = [{k: v for k, v in r.items() if k != "pediatric_context"} for r in rows]
    assert json.dumps(stripped, sort_keys=True) == before, "annotation changed an existing field"
    doc["pediatric_context"] = {
        "note": annotate.NOTE,
        "provenance": inp.provenance(),
        "min_samples_per_group": min_n,
        "lof_classes": sorted(annotate.LOF_CLASSES),
        "dna_samples": "one WGS/WXS primary tumour per patient (targeted panels excluded)",
        "genes": [{k: v for k, v in rec.items() if not (isinstance(v, float) and math.isnan(v))}
                  for rec in ctx.to_dict(orient="records")],
    }
    Path(out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(out_json).write_text(json.dumps(doc, indent=1, allow_nan=False) + "\n")
    log(f"annotated {len(rows)} rows -> {out_json}")
    return doc
