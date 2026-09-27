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

## Follow-up run: health-referenced patients, held-out signature, retrieval check

The four next steps from the first run are done.

**1. Healthy reference without cross-batch z-scores.** GSE244118 (abdominal subcutaneous
adipose RNA-seq, Klein lab) contrasts obese (MHO + MUO) with metabolically healthy lean
(MHL) adults within one batch. `build-offset` turns that contrast into a per-gene
Cohen's d. The 4 `PE`-prefix libraries were all MUO, so the offset uses the `PSQ`
batch only: 34 obese vs 15 lean. The patient's z against the obese GSE95640 cohort,
plus this offset, is their z against health:
`(x − μ_lean)/σ = (x − μ_obese)/σ + (μ_obese − μ_lean)/σ`. This assumes the
between-subject SD is similar in both studies.

The offset recovers known obesity markers:

| direction | genes (Cohen's d) |
|---|---|
| up in obese | LEP +2.4, SPP1 +2.5, CD68 +2.4, MMP9 +2.3, SAA1 +1.8 |
| down in obese | CIDEA −3.0, ADH1B −2.9, GPX3 −1.3 |

On the top 500 obesity genes, the LCD signature anti-correlates with it (cos −0.42).
Two independent cohorts agree that weight loss moves adipose back toward lean.

**2. Held-out signature.** `ingest-npi --exclude-samples` drops the patient's whole
subject before building the signature.

**3. Drug-side tissue penalty.** A new `drug_tissue_mismatch` shrink (×0.8) applies to
every LINCS cell-line signature scored against a tissue. It is uniform across LINCS
drugs, so it lowers confidence without changing their order.

**4. Retrieval regression check** (`scripts/retrieval_check.py`). The paper's
L1000CDS2 top 50 mimics come from Supplementary Table S5 (Phase I signatures; now at
`configs/retrieval/GSE95640_L1000CDS2_top50.tsv`). 18 of its 43 compounds are in
Phase II. Each of the 1,796 Phase II compounds was ranked by cosine to our LCD
signature:

| per-compound aggregate | AUC of the paper's hits | permutation p | paper hits in our top 10% |
|---|---|---|---|
| median over signatures | 0.62 | 0.04 | 7/18 |
| best signature (L1000CDS2-like) | 0.75 | 0.0005 | 7/18 |

The replicated hits are the PI3K/mTOR inhibitors: GSK-1059615, NVP-BEZ235,
WYE-125132 and GSK-2126458, all in our top 2%. Tanespimycin and MG-132 also
replicate. Found independently, our top 15 mimics are all PI3K/mTOR inhibitors. That
fits LCD biology: less insulin/mTORC1 → SREBP1 signalling, and lower SCD/FADS1/FADS2.
The paper's HDAC-inhibitor hits do not replicate (trichostatin A is in the 1st
percentile, vorinostat in the 45th).

### Cohort sweep: all 191 subjects as patients (`scripts/cohort_sweep.py`)

Each subject was held out of the LCD signature, encoded against health via the
offset, and scored against the 31 drugs.

| quantity | result |
|---|---|
| LCD reverses s_P (reverse_npi > 0) | 79.6% of subjects (median +0.34) |
| subjects with any drug of positive complementarity | 68.6% |
| metabolic set vs random background, per-patient AUC | median 0.495 (IQR 0.39–0.57) |
| between-patient Spearman of drug ranks | median 0.23 |
| share of \|s_P\|² carried by the shared offset | median 0.49 |
| most frequent top named drug | everolimus (67), berberine (48), resveratrol (24), metformin (19) |

What this shows:

- **Patient side:** it now behaves as intended. LCD reverses most subjects'
  health-referenced state. Rankings differ between patients rather than just
  repeating the shared offset.
- **Drug side:** it is the bottleneck. The metabolic drugs do not outrank random
  compounds, and random background compounds (VT-464, aripiprazole) reach rank 1 as
  often as rapalogs. L1000 cancer-cell-line consensus signatures do not carry adipose
  pharmacology well enough to tell these drugs apart.
- **Where the drug data does agree:** the one class with consistent signal is
  mTOR/PI3K inhibition. It is the top retrieved LCD mimic, and everolimus is the most
  frequent top named drug. So "which drug adds to LCD" and "which drug is LCD-like"
  currently land on the same class. That points to redundancy, not complementarity,
  and scores should be read with that in mind.

### Next

- Use adipose- or adipocyte-derived drug signatures where they exist (e.g. the Phase I
  ASC and HA1E plates behind the paper's hits, or adipocyte perturbation datasets)
  instead of pan-cancer consensus.
- Get outcome data, such as weight regain, to test whether patient-level scores
  predict anything. GSE95640 withholds clinical metadata.

## Reproduce

Follow-up run:

```bash
npi-pharma fetch-geo GSE244118 --suppl && python scripts/curate_gse244118.py
npi-pharma build-offset --counts data/raw/GSE244118/GSE244118_abdominal.fat_all.gene_counts.txt.gz \
  --id-map data/raw/ensembl_to_symbol.tsv --samples data/raw/GSE244118/samples_curated.tsv \
  --filter batch=PSQ --case MHO,MUO --control MHL --out data/processed/offset_obese_vs_lean_GSE244118.tsv
npi-pharma ingest-npi --ids LCD_adipose_GSE95640 --exclude-samples B00EVC7 --out npis_holdout.parquet
npi-pharma encode-patient --input baseline_logcpm.tsv --sample B00EVC7 --tissue adipose \
  --reference-offset data/processed/offset_obese_vs_lean_GSE244118.tsv --out patient.npz
python scripts/cohort_sweep.py --offset data/processed/offset_obese_vs_lean_GSE244118.tsv   # ~2 min
python scripts/retrieval_check.py                                                         # ~4 min, streams the GCTX
```

First run:

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
