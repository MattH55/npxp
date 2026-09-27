# First real-data milestone: GSE95640 LCD × LINCS L1000

Run on 2026-09-27. Every input is real public data; nothing is synthetic. All
scores are `tier_3_mechanism_only`.

## Inputs

| input | source | notes |
|---|---|---|
| NPI | GEO GSE95640 (DiOGenes), subcutaneous adipose RNA-seq, 191 obese non-diabetic adults, baseline vs after an 8-week LCD at 800–1000 kcal/d | Raw counts come from the supplementary file. They are mapped Ensembl → symbol via NCBI gene_info and converted to log2(CPM+1), keeping genes with CPM ≥ 1 in ≥ 50% of samples (14,396 genes). The contrast is paired post − pre, n = 191. |
| drugs | GEO GSE70138 (LINCS Phase II) Level 5 COMPZ, per-compound median across all cell lines and doses | 11 of the 12 `metabolic` drugs were found (fenofibrate is absent from Phase II and reported as missing). A background of 20 random `trt_cp` compounds was added with seed 0. |
| patient | GSE95640 baseline sample B00EVC7 (subject S001) | There is no healthy-adipose reference on this platform. s_P is taken against the other 190 baselines and flagged `cohort_relative`. S001 is also 1 of the 191 pairs behind the NPI signature. |

### GSE95640 curation

GEO gives no subject id. `scripts/curate_gse95640.py` pairs adjacent rows of the
supplementary sample sheet and refuses to write the curated sheet unless three
checks pass:

- **Pairing, metadata:** all 191 pairs have one CID1 and one CID2 sample and matching sex and age.
- **Pairing, expression:** the median within-pair expression correlation is 0.46, against 0.16 for shifted pairs.
- **Timepoint direction:** the series matrix says "CID1 = baseline, CID2 = after 8 weeks of LCD". The supplementary sheet's header says "CID1 = after LCD, CID2 = 6 months after". The data agree with the series matrix: SCD falls by −1.66 log2 from CID1 to CID2.

### NPI signature sanity check

The signature matches published LCD adipose biology:

- **Top down:** SCD (z −9.3), FADS2, FADS1, ALDOC, CES1 (de novo lipogenesis and desaturation).
- **Top up:** C6, CIDEA, GPX3, ADH1B, SLC27A2 (genes reported to rise after weight loss).

## Result (patient B00EVC7, `rank --npi-class diet --drug-set metabolic`)

| rank | drug | composite | complementarity | rev(drug) |
|---|---|---|---|---|
| 1 | berberine | −0.071 | −0.170 | +0.41 |
| 2 | resveratrol | −0.073 | −0.176 | +0.38 |
| 3 | acarbose | −0.089 | −0.153 | |
| 5 | sirolimus | −0.096 | −0.095 | |
| 8 | metformin | −0.123 | −0.184 | −0.30 |
| 11 | simvastatin | −0.186 | −0.258 | |

## What this does and does not show

- **LCD does not reverse this subject's state.** rev(LCD) is −0.22: this subject's
  deviation from the obese cohort already points the same way as the LCD response.
  Complementarity is therefore negative for every drug, and the ranking is driven
  by the drug's own reversal of s_P. In this configuration the table is a ranking
  of single drugs, not of NPI × drug complementarity.
- **The metabolic set is not enriched over the random background.** Scored against
  LCD, the AUC for a metabolic drug outranking a random compound is 0.32. On this
  input the scores separate nothing. That is the expected outcome for a
  cohort-relative s_P with no disease contrast.
- **Drug-side tissue mismatch is not penalised.** L1000 signatures come from cancer
  cell lines, while the patient and NPI are adipose. `confidence` stays ×1.00 because
  the tissue shrink factor is applied only to NPI vs patient tissue. This is a gap in
  the scoring, not a property of the data.

### Next steps for a meaningful per-patient test

1. Use a healthy (lean) adipose reference on a comparable RNA-seq pipeline, such as
   GTEx subcutaneous adipose, so that s_P is a disease vector rather than a cohort
   deviation.
2. Hold the patient's own pair out of the NPI signature.
3. Add a drug-side cell-line-vs-tissue shrink factor.
4. Curate the L1000CDS2 hits reported in the GSE95640 paper and use them for the
   `retrieval_references` regression check.

## Reproduce

```bash
npi-pharma fetch-geo GSE95640 --suppl --gene-info
python scripts/curate_gse95640.py
npi-pharma ingest-npi --ids LCD_adipose_GSE95640
npi-pharma fetch-lincs --gse GSE70138        # 5.4 GB download, 5.8 GB GCTX on disk
npi-pharma ingest-lincs --gctx ... --drug-set metabolic --random 20 --seed 0 --out data/processed/signatures/drugs.parquet
# patient: GSE95640 CID1 samples as a genes x samples log-CPM TSV
npi-pharma encode-patient --input baseline_logcpm.tsv --sample B00EVC7 --tissue adipose --out patient.npz
npi-pharma rank --patient patient.npz --npi-class diet --drug-set metabolic --top 15 --out rank.tsv --json rank.json
```
