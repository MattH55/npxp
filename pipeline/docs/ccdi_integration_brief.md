# Coding-agent brief: integrate CCDI pediatric data sources into NPxP

Source list: the NCI CCDI slide "Combined Data Drives Scientific Discovery"
(Childhood & AYA Data Innovation Challenge, 2026). This brief says which of
those sources to integrate, in what order, and into which part of NPxP. It
also says what must not be built.

> **Status (2026-10-20):** milestones 1–3 have been run on OpenPedCan v15.
> Results and deviations from this brief are in `ccdi_integration.md`.
> Notably, §4.2's joint-PC1 concordance recipe failed its positive control
> and was replaced.

Work in `pipeline/` (`npi_pharma`). Read these first:

- `pipeline/README.md`
- `pipeline/docs/milestone_GSE95640.md`
- `pipeline/docs/validation_drugcomb.md`
- `pipeline/docs/mechanism_hypotheses.md`
- `BUILD_SPEC.md` (the `cell_line` schema, especially `oncotree_*` and
  `tumor_concordance_score`)

---

## 0. Ground rules (non-negotiable)

1. **Do not build a pediatric synergy predictor.** `validation_drugcomb.md`
   showed that signature composition has no measurable information about
   synergy on 739,964 real drug pairs. CCDI data must not be used to retrain,
   re-weight or re-advertise `interact/score.py` as predictive. What it is used
   for:
   - patient and tumour **context**
   - **model-fidelity** checks
   - **annotations** on the existing mechanism hypotheses

   Every new output carries an evidence tier, and none can read as more
   confident than `tier_3_mechanism_only`.
2. **Open and registered access only, in code.**
   - **Open:** GDC open tier (TARGET, HCMI), refine.bio, Treehouse public
     compendium, OpenPedCan / Molecular Targets Platform, ScPCA, FDA PMTL,
     SEER*Explorer statistics.
   - **Login token:** PedcBioPortal gets an optional adapter that reads
     `PEDCBIO_TOKEN` from the environment.
   - **Controlled (dbGaP / DUA):** Kids First, St. Jude Cloud, COG / NCTN,
     Project:EveryChild, CCSS, Foundation Medicine, PPCR. These get a
     documentation entry only (§7): no code, no data, no derived per-patient
     output.
3. **No fabricated numbers.** If a file, field or endpoint is not what this
   brief says, stop and record what you actually found. Do not substitute a
   guess. Every API detail below is marked **(verify)** because it was written
   without network access. Confirm each one against the live service before
   coding against it.
4. **Provenance on every artefact.** Each download writes
   `data/raw/<source>/MANIFEST.json` with:
   - source, URL, release/version and retrieval date (UTC)
   - sha256 and size of each file
   - access tier and licence / citation text

   Every derived table carries the manifest's release and version.
5. **Never commit raw data or per-patient output.** `data/` and `out/` are
   git-ignored. The exceptions are small curated tables under
   `pipeline/configs/ccdi/` (e.g. the PMTL list, SEER summary statistics), each
   with a source URL and retrieval date in a header comment.
6. **`api/` is generated.** Files exported by the backend's `export_static.py`
   (`opportunity_*.json`, `cell_lines*`, `drugs*`, `interactions*`, …) must not
   be hand-edited. Pipeline-produced pages follow the
   `api/mechanism_hypotheses.json` precedent: the pipeline script writes to
   `out/`, and the JSON is copied into `api/` together with its page. Anything
   that changes a backend table goes into a hand-off note (§6) instead.
7. **No network in unit tests.** Use tiny fixtures in the real file formats
   (the same pattern as `tests/test_fetch.py` and `tests/test_ingest.py`).
   Network-dependent steps belong in `scripts/` or behind CLI commands.
8. **Personal data.** refine.bio requires an email address for dataset
   requests. Read it from `REFINEBIO_EMAIL`. Never hard-code it, and never
   default it to any address found in git config.

---

## 1. Network and access preflight (do first)

Add `scripts/ccdi_preflight.py`. It probes each host below and prints a
reachable / blocked table:

- `api.gdc.cancer.gov`
- `api.refine.bio`
- `treehousegenomics.soe.ucsc.edu`
- the OpenPedCan S3 bucket host
- `moleculartargets.ccdi.cancer.gov`
- `api.scpca.alexslemonade.org`
- `pedcbioportal.kidsfirstdrc.org`
- `www.fda.gov`
- `seer.cancer.gov`

If a host is blocked, report it and continue with the sources that are
reachable. Do not look for mirrors of unknown provenance.

---

## 2. New package layout

```
pipeline/src/npi_pharma/ccdi/
  __init__.py
  manifest.py      # MANIFEST.json read/write, sha256, retrieval date
  gdc.py           # TARGET + HCMI open-tier STAR counts
  refinebio.py     # dataset request/poll/download (token + REFINEBIO_EMAIL)
  treehouse.py     # public compendium + clinical table
  openpedcan.py    # histologies, mutation/CNV/fusion frequencies, TPM (or MTP API)
  scpca.py         # project listing + h5ad download (token)
  pedcbio.py       # optional, PEDCBIO_TOKEN
  oncotree.py      # diagnosis string -> OncoTree code, one table for every source
pipeline/configs/ccdi/
  fda_pmtl.tsv     # curated, dated
  seer_pediatric.tsv
  diagnosis_map.tsv  # source diagnosis label -> OncoTree code (curated, reviewed)
```

**CLI.** Add one command per source: `fetch-gdc`, `fetch-refinebio`,
`fetch-treehouse`, `fetch-openpedcan` and `fetch-scpca`. Each writes to
`data/raw/<source>/` plus a MANIFEST. Add one `build-ccdi` command that
derives the tables in §4.

**Harmonisation, one function each, shared by every source:**
- **Genes:** map to HGNC symbols with `genes.normalize_symbol` /
  `dedupe_vector`. For Ensembl IDs, strip the version and map through the
  existing NCBI gene_info map (`fetch.fetch_gene_info`).
- **Expression:** convert to log2(TPM+1). Where only counts exist, use
  log2(CPM+1) and say so in the MANIFEST.
- **Diagnosis:** map to OncoTree through `configs/ccdi/diagnosis_map.tsv`.
  Unmapped labels are reported, never dropped silently.
- **Age:** store `age_at_diagnosis_years` and `pediatric_aya` (< 40, with < 18
  and 18–39 kept distinct). Never infer age from the cohort name.

---

## 3. Sources: what to fetch (all details **(verify)**)

### 3.1 GDC: TARGET and HCMI (open tier)
- **File search:** `POST https://api.gdc.cancer.gov/files`.
  - Filters:
    - `cases.project.project_id` in the TARGET projects. List them live: do not
      hard-code. Expected examples are TARGET-NBL, TARGET-WT, TARGET-OS,
      TARGET-AML and TARGET-ALL-P2, plus `HCMI-CMDC`.
    - `data_type = "Gene Expression Quantification"`
    - `analysis.workflow_type = "STAR - Counts"`
    - `access = "open"`
  - Fields: `file_id`, `file_name`, `md5sum`, `cases.submitter_id`,
    `cases.samples.sample_type`, `cases.project.project_id`,
    `cases.diagnoses.age_at_diagnosis` (in days),
    `cases.diagnoses.primary_diagnosis`.
- **Download:** `POST https://api.gdc.cancer.gov/data` with up to N ids per
  request, which returns a tar.gz. Verify each file's md5 against the search
  result.
- **File format:** TSV with a `# gene-model:` comment line, then
  `N_unmapped` / `N_multimapping` / … summary rows to drop. Columns include
  `gene_id`, `gene_name`, `gene_type`, `unstranded` and `tpm_unstranded`. Keep
  `protein_coding` plus lncRNA, and record which.
- **Sample types:**
  - Keep `Primary Tumor`, `Recurrent Tumor` and the normal types
    (`Solid Tissue Normal`, `Blood Derived Normal`, `Bone Marrow Normal`)
    separate.
  - HCMI models carry model sample types such as "Next Generation Cancer
    Model". Tag them `kind=model`.
  - Pediatric HCMI is the subset with age_at_diagnosis < 21 years. Report the
    count you find; do not assume one.

### 3.2 refine.bio (Childhood Cancer Data Lab)
- **Token:** `POST https://api.refine.bio/v1/token/`, then activate it (accept
  the terms) by `PUT`-ing `is_activated: true`.
- **Dataset request:** `POST /v1/dataset/` with
  `{"data": {"<accession>": ["ALL"]}, "email_address": $REFINEBIO_EMAIL, "aggregate_by": "EXPERIMENT", "scale_by": "NONE", "quantile_normalize": false}`.
  Then `PUT /v1/dataset/{id}/` with `{"start": true}`, poll until
  `is_processed`, and download `download_url`.
- **Use:**
  - Uniformly processed versions of GEO series already in
    `configs/npi_catalog.yaml`, so they can be cross-checked against the
    hand-parsed series.
  - Pediatric tumour microarray cohorts that have no RNA-seq.
- **Do not mix values across platforms in one z-score.** refine.bio
  microarray and RNA-seq values are on different scales.

### 3.3 Treehouse public compendium
- **Files:** discover them from the public-data page; do not guess file names.
  Expect a log2(TPM+1), HUGO-symbol compendium plus a clinical table
  (disease, age, pediatric / AYA flag).
- **Use:**
  - The cross-cohort pediatric tumour **cohort reference** for
    `encode-patient`, per disease.
  - The concordance target in §4.2.

### 3.4 OpenPedCan / Molecular Targets Platform
- **Data.** Prefer the Molecular Targets Platform (MTP, Open Targets-based)
  API for precomputed per-histology aggregates, if it exists and is reachable:
  - SNV / CNV / fusion frequencies per gene × histology
  - TPM summaries
  - PMTL designation

  Otherwise use the OpenPedCan release files from the public bucket, as listed
  by the OpenPedCan-analysis repo's download script: `histologies.tsv`, the
  consensus SNV MAF, consensus CNV, putative-oncogenic fusions and the TPM
  matrix.
- **The TPM matrix may be an `.rds`.** Use `pyreadr` if it can read it.
  Otherwise follow the repo's own precedent for unreadable RDS files (the R
  workaround in the CTRP commit) and document it.
- **Use:** pediatric context annotations on the mechanism hypotheses (§4.3).

### 3.5 ScPCA portal (single-cell)
- **Access:** list projects via the API and download per-sample AnnData
  (`.h5ad`) with a token.
- **Use:** feed `patient/encode.read_expression(..., sample_key=..., obs_filter=...)`.
  Pseudobulk tumour cells only, using ScPCA's own cell-type and tumour
  annotations. Record which annotation column was used.

### 3.6 FDA Pediatric Molecular Target Lists (RACE for Children Act)
- **Data:** a curated TSV, `configs/ccdi/fda_pmtl.tsv`, with columns gene,
  list (relevant / non-relevant), FDA list version and retrieval date. If MTP
  exposes PMTL, cross-check the two and report the differences.
- **Use:** an annotation on drugs and targets. It does not change scores.

### 3.7 SEER (and NPCR) pediatric statistics
- **Data:** `configs/ccdi/seer_pediatric.tsv`, a SEER*Explorer export of
  incidence per 100k and 5-year relative survival for ages 0–19, by ICCC site.
  Map it to OncoTree through `diagnosis_map.tsv`. Cite it as SEER requires.
- **Use:** hand-off only (§6), for weighting the Opportunity Ranking by burden.

### 3.8 PedcBioPortal (optional)
- **Data:** cBioPortal REST under `/api`, with a bearer token from
  `PEDCBIO_TOKEN`. Skip the source entirely when the token is unset.
- **Use:** a fallback for studies that OpenPedCan / MTP does not aggregate.

---

## 4. Derived outputs (`npi-pharma build-ccdi`)

### 4.1 Pediatric patient context (`encode-patient`)
- **New mode:** `--cohort-ref data/processed/ccdi/<disease>.parquet`.
- **s_P definition:** the patient's z-score against tumours of the **same
  OncoTree code**, flagged `cohort_relative`. Real pediatric tumours have no
  healthy-tissue reference. Do not z-score tumours against GDC normals from
  another tissue.
- **Tissues:** add pediatric canonical tissues to `scoring.yaml` (e.g.
  `peripheral_nervous_system`, `bone`, `kidney`, `bone_marrow`, `cns`), so
  that `tissue_mismatch` against the adipose / muscle NPI catalog fires
  correctly. Expect it to fire for almost every current NPI. That is correct
  behaviour: report it, and do not loosen the synonyms to avoid it.

### 4.2 Model fidelity: fill `tumor_concordance_score`
For every NPxP cell line, and every pediatric HCMI or DepMap model, against
the tumour cohort of its OncoTree code (TARGET + Treehouse + OpenPedCan):

1. Take the top 1,000 most variable genes within that tumour cohort.
2. Remove the first principal component from tumours and models jointly, which
   is a crude purity / batch adjustment. Document this as an approximation of
   Celligner, not a replication of it.
3. `concordance` = Spearman correlation of the model's profile with the cohort
   centroid.
4. Report a percentile against a null of models from *other* OncoTree codes.
5. Output `out/ccdi/model_concordance.tsv` with columns `cell_line_id`,
   `oncotree_code`, `n_tumours`, `concordance`, `null_percentile`, `flag`
   (`concordant`, `discordant`, `no_cohort`).
6. Lines with no pediatric cohort get `no_cohort`, not a number.

### 4.3 Pediatric annotations on mechanism hypotheses
For each row of `out/mechanism_hypotheses.json`, add a `pediatric_context`
block without changing any existing field. For each OncoTree code with
≥ 20 tumours, it reports:

- **Baseline lesion frequency.** For each lesion or SL-partner gene the row
  relies on (e.g. BRCA1/2 and other HR genes for the induced-HRD rows), the
  fraction of tumours with a putative loss-of-function SNV, deep deletion or
  fusion. This says whether the tumour already carries the state the NPI is
  supposed to induce.
- **Target expression.** Median and IQR log2(TPM+1) of each drug's target
  (e.g. PARP1) per histology, plus its percentile across all pediatric
  histologies.
- **Regulatory flag:** the PMTL designation of the target.
- **Pediatric trials:** reuse `scripts/fetch_clinical_trials.py` with an age
  filter (minimum age < 18). List NCT IDs only.

These are lookups, each with its source release. Add no composite or ranking
on top of them. Re-export `api/mechanism_hypotheses.json` and update
`js/mechanism.js` to show the block as a collapsible section per row,
labelled "Pediatric context (annotation, not evidence of interaction)".

### 4.4 Cross-check of refine.bio against hand-parsed GEO
For each `npi_catalog.yaml` series that refine.bio also has, rebuild the NPI
signature from refine.bio values and report the cosine against the existing
signature on shared genes. This is a regression check on the curation scripts:
a cosine below 0.8 is investigated, not ignored.

---

## 5. Tests (offline; tiny fixtures in real formats)

- **GDC:**
  - The search payload builder emits the expected filters.
  - The STAR-counts parser drops the `N_*` rows and the comment line, and maps
    versioned Ensembl IDs.
  - An md5 mismatch raises.
- **refine.bio:**
  - Request / poll / download against a mocked `urlopen`.
  - The client refuses to run when `REFINEBIO_EMAIL` is unset.
- **Treehouse / OpenPedCan:** a histology filter at the ≥ 20 tumours
  threshold. LoF frequency on a 6-row MAF fixture, including a non-LoF variant
  that must not count.
- **OncoTree mapping:** unmapped labels are returned in a report, not dropped.
- **Concordance:**
  - A synthetic model built from its own cohort centroid plus noise scores
    `concordant`.
  - A model from another cohort scores below the 50th null percentile.
  - No cohort gives `no_cohort`.
- **Mechanism annotation:** existing fields are byte-identical before and after
  the annotation is added.
- **Manifest:** sha256 and retrieval date are written. Re-running without
  `--overwrite` does not re-download.
- **Patient context:** a tumour patient with an adipose NPI gets
  `tissue_mismatch`.

---

## 6. Hand-off note for the backend (not in this repo)

Write `pipeline/docs/ccdi_backend_handoff.md`. It lists the proposed changes
to backend-exported tables, with the file that supplies each value:

- `cell_line.tumor_concordance_score` and `tumor_concordance_flag`, from
  `model_concordance.tsv`.
- New pediatric `cell_line` rows (HCMI / DepMap pediatric models), with
  OncoTree codes.
- Opportunity Ranking: optional burden columns from `seer_pediatric.tsv`.
  Leave the existing score formula unchanged unless a separately documented
  methodology change is approved.

Do not edit `api/opportunity_*.json` or `api/cell_lines*` by hand.

---

## 7. Controlled-access sources: documentation only

In `pipeline/docs/references.md`, add a "CCDI controlled-access" section. For
each source below, give one line covering the data type, the access route
(dbGaP study / DUA / application) and what NPxP could use it for later. Write
no code.

- Kids First DRC
- St. Jude Cloud
- COG / Pediatric MATCH / MyPART (NCTN Data Archive)
- Project:EveryChild
- CCSS. Note it as the one source with lifestyle (exercise / diet) × treatment
  × late-effect data. It would need a new population-level evidence tier,
  outside the current cell-line schema.
- Foundation Medicine
- PPCR

---

## 8. Milestones and acceptance

1. **Preflight + manifests + GDC TARGET-NBL only.** End to end: fetch,
   harmonise, then `encode-patient --cohort-ref` for one neuroblastoma sample.
   Done when the tests pass and the MANIFEST is complete.
2. **Concordance** for every NPxP cell line that has a pediatric OncoTree
   cohort, plus pediatric HCMI. Done when the TSV and a short results section
   are written. Results must report the discordant lines too, not only the
   good ones.
3. **Mechanism-hypothesis pediatric annotations**, the PMTL table and
   pediatric trials. Done when the page re-export and `js/mechanism.js` change
   are in, with existing hypothesis fields unchanged.
4. **refine.bio cross-check (§4.4)**, the ScPCA pseudobulk path and the
   backend hand-off note.

After each milestone, append real numbers, including null or negative
results, to `pipeline/docs/ccdi_integration.md`. Put the commands, releases and
retrieval dates in `pipeline/README.md`. Stop after milestone 1 if a source in
§1 is blocked, and report which.
