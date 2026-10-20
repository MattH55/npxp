# CCDI integration: results

Run on 2026-10-20 in a cloud session, following `ccdi_integration_brief.md`. All
inputs are real, open data. All outputs are annotations or model-fidelity
checks, not interaction evidence.

## 1. Which sources could be reached

`scripts/ccdi_preflight.py`:

| source | status |
|---|---|
| **OpenPedCan public bucket** (`s3.amazonaws.com/d3b-openaccess-us-east-1-prd-pbta`) | **reachable**, and used for everything below |
| GDC, refine.bio, Treehouse, Molecular Targets Platform, ScPCA, PedcBioPortal, FDA, SEER, NCBI GEO | blocked by the session's network policy (proxy 403) |

**OpenPedCan stands in for most of the blocked sources.** It re-processes
several cohorts through one shared pipeline:
- TARGET (the TARGET data otherwise on GDC)
- PBTA/CBTN
- GMKF (Kids First)
- PPTC xenografts
- the Maris lab neuroblastoma cell lines

It also ships the FDA Pediatric Molecular Target List mapping
(`ensg-hugo-pmtl-mapping.tsv`, PMTL v3.1), which stands in for the blocked FDA
page.

**Still missing:** Treehouse, refine.bio, ScPCA, SEER and HCMI. Those steps
(§4.4, single-cell, burden weighting) have not been run.

### Release and files

OpenPedCan **v15**, the last complete release; v16 is an incremental drop with
no `histologies.tsv`. Licence CC-BY-4.0. Every file is md5-verified against the
release's `md5sum.txt`, and sha256, size and retrieval time are recorded in
`data/raw/openpedcan/v15/MANIFEST.json`.

| file | MB |
|---|---|
| `histologies.tsv` | 21 |
| `ensg-hugo-pmtl-mapping.tsv` | 2 |
| `independent-specimens.rnaseqpanel.primary.tsv` | 1 |
| `independent-specimens.wgswxspanel.primary.prefer.wgs.tsv` | 0.3 |
| `gene-expression-rsem-tpm-collapsed.rds` | 285 |
| `fusion-putative-oncogenic.tsv` | 18 |
| `consensus_wgs_plus_cnvkit_wxs_plus_freec_tumor_only.tsv.gz` | 510 |
| `snv-consensus-plus-hotspots.maf.tsv.gz` | 2,077 |

The TPM `.rds` is an R data.frame of 60,325 genes × 4,124 samples. `pyreadr`
reads it in about 50 s and 4 GB of RAM, so R is not needed. It is converted
once to a log2(TPM+1) parquet.

```bash
npi-pharma fetch-openpedcan --with-maf
npi-pharma build-ccdi                      # cohorts, concordance, annotate
```

## 2. Pediatric tumour cohorts

The cohorts are built from one primary RNA specimen per patient, using
OpenPedCan's independent-specimen list. TCGA (adult) samples are excluded.

**Cohorts are stratified by RNA library, not just by cancer group.** The
neuroblastoma group, for example, mixes TARGET poly-A libraries (123 tumours)
with GMKF stranded libraries (167). A z-score against that mix mostly measures
library chemistry. The result is 31 cohorts (cancer group × library) with
≥ 20 tumours, listed in `docs/ccdi/cohort_strata_openpedcan_v15.tsv`.

**OncoTree mapping.** `configs/ccdi/diagnosis_map.tsv` is hand-curated, and
its codes were not checked live because oncotree.info was blocked. 45 cancer
groups are left unmapped and reported. Some are deliberate: broad groups such
as "Low-grade glioma" span several codes.

### Milestone 1: one neuroblastoma patient

TARGET-30-PAIFXV-01A-01R (MYCN non-amplified) was encoded against the TARGET
poly-A neuroblastoma cohort, n = 122 after leaving the patient out. It is
flagged `cohort_relative`: s_P is the deviation from other neuroblastomas,
not from health. Real pediatric tumours have no healthy-tissue reference, and
GTEx adrenal tissue (mostly cortex, a different batch) is not an acceptable
substitute.

| check | result |
|---|---|
| MYCN z (single patient) | −1.71 |
| E2F proliferation program (single patient) | −11.7 |
| inflammatory program (single patient) | +3.8 |
| NTRK1 z (single patient) | −1.4, unexpected: high NTRK1 usually goes with MYCN non-amplified tumours |
| MYCN z across the cohort, amplified (n = 63) vs non-amplified (n = 60) | median +0.18 vs −0.50, AUC 0.81 |

The cohort-level AUC shows a real separation, weaker than expected for an
amplification marker. Scored against the adult adipose / muscle NPI catalog,
this patient gets `tissue_mismatch`, which is the correct behaviour. It is
covered by a test.

## 3. Model fidelity: models vs pediatric tumours

The output is `docs/ccdi/model_concordance_openpedcan_v15.tsv`: 489 models,
comprising 206 cell lines and 283 PDXs.

### Final method

1. Run the comparison separately per model type (cell line, PDX).
2. The reference groups are the cancer groups with ≥ 20 tumours *and* ≥ 1
   model of that type.
3. Genes: the 1,000 expressed genes with the highest variance between those
   groups' tumour centroids.
4. Centre tumour centroids on their mean, and models on the mean of their
   per-group means.
5. Concordance is the Spearman correlation of a model with its own group's
   centroid. The null is same-type models from other groups.
6. `concordant` means the model beats ≥ 95% of the null *and* its own group is
   its best match.

| | concordant | intermediate | discordant | no tumour cohort |
|---|---|---|---|---|
| cell lines (10 reference groups) | 51 | 48 | 13 | 94 (87 chordoma) |
| PDXs (9 reference groups) | 194 | 11 | 1 | 77 |

| models | concordant | median ρ |
|---|---|---|
| neuroblastoma cell lines (Maris) | 33 / 38 | 0.57 |
| neuroblastoma PDXs | 33 / 33 | 0.67 |
| osteosarcoma PDXs | 30 / 30 | 0.59 |
| Wilms PDXs | 11 / 12 | 0.60 |
| medulloblastoma PDXs | 15 / 18 | 0.52 |
| ependymoma PDXs | 6 / 6 | 0.56 |
| **high-grade glioma cell lines** | **5 / 32** (8 discordant) | **−0.14** |
| diffuse midline glioma cell lines | 4 / 20 | 0.22 |

- **Glioma cell lines are the weakest models.** That is consistent with the
  known drift of glioma lines in serum culture.
- **The ALL PDX result (90/90) is close to trivial.** ALL is the only
  leukaemia among the PDX reference groups, so "matches ALL best" partly means
  "is a leukaemia".

### How the method was arrived at

It took three versions, and the brief's own §4.2 recipe failed. The Maris
neuroblastoma lines, which are established neuroblastoma models, served as the
positive control throughout.

1. **Brief's recipe:** remove PC1 jointly from tumours and models. Every
   neuroblastoma line came out *anti*-correlated with neuroblastoma tumours
   (ρ ≈ −0.5), with AML as the best match. Joint centring makes the tumour
   centroid mostly encode tumour vs culture.
2. **Separate centring of tumours and models:** xenografts were fixed, but the
   neuroblastoma lines still matched ALL best. Tumour centroids were centred
   over groups that include the leukaemias, while the cell lines are all solid
   tumours.
3. **Per model type over shared groups:** the final method above.

Because the method was revised against a known-positive control, the flags
should be read as a screen. `tests/test_ccdi.py` checks the final method on a
synthetic world with a large shared culture offset, a mislabelled model and a
model with no cohort.

**NPxP's own 16 cell lines are not assessable here.** None is a pediatric
model, and none has expression in OpenPedCan. DepMap was not reachable. Their
`tumor_concordance_score` stays null.

## 4. Pediatric context on the mechanism hypotheses

`api/mechanism_hypotheses.json` was re-exported with:
- a `pediatric_context` block on each of the 614 rows;
- a top-level `pediatric_context.genes` table (68 genes × 22 cancer groups,
  1,496 records).

Every existing field is byte-identical to the previous export, which is
asserted in code. `mechanism-hypotheses.html` shows the blocks as collapsible
sections labelled "Pediatric context (annotation, not evidence of
interaction)".

**DNA denominators** are one WGS/WXS primary tumour per patient. Targeted
panels are excluded, so a gene outside a panel cannot read as "never
mutated". LoF means nonsense, frameshift, splice-site, nonstop or
start-loss; missense variants are not counted. Somatic calls only: germline
predisposition is not included.

**Baseline HR lesions are rare in every pediatric group.** For the 590
induced-HR-deficiency rows, no group reaches ≥ 5% baseline LoF or deep
deletion in the lesion genes:

| gene | highest group | altered fraction |
|---|---|---|
| BRCA2 | high-grade glioma | 3.6% (4 / 112) |
| RAD51 | ganglioneuroblastoma | 1.4% |
| BRCA1 | ALL | 0.4% |

So the hyperthermia / hypoxia / TTFields "induced BRCAness" hypotheses would
act on tumours that are almost all HR-proficient at baseline. That is the
premise of those hypotheses, now measured in pediatric cohorts.

**Fasting/CR and exercise rows** (24) get no lesion frequency: their "lesion"
genes are pathway or immune targets (PIK3CA, MTOR, PDCD1, …), where LoF
frequency is meaningless. They get target expression and PMTL only.

**Targets:**
- PARP1 is most expressed in Wilms tumour (median log2 TPM 7.3), then diffuse
  midline glioma and ATRT.
- PMTL: 44 of the 68 genes are FDA "Relevant Molecular Targets", 1 is
  "Non-Relevant", and 23 are not listed.

**Not done:** pediatric trial lists, because ClinicalTrials.gov was not
re-queried in this session.

## 5. What remains

| item | blocked on |
|---|---|
| GDC-native TARGET / HCMI, Treehouse, refine.bio cross-check (§4.4) | network access to those hosts |
| ScPCA single-cell pseudobulk | network + token |
| SEER burden columns, FDA PMTL cross-check against the live list | network |
| concordance for NPxP's 16 adult lines | DepMap expression (backend workspace) |
| backend table changes | `ccdi_backend_handoff.md` |
