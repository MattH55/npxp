# NPI × drug predictor for cancer cell lines: resistance-program reversal

Run 2026-09-28. Code: `src/npi_pharma/cancer/`, `scripts/cancer_npi_drug.py`.

**Verdict: the predictor does not reproduce the measured NPI × drug results, so its
ranking is not usable yet.** The drug side works; the step from NPI signature
to "sensitises to drug" does not.

## Method

- **Drug resistance programs.** For each drug and each of the 949 LINCS landmark
  genes, the correlation across cancer cell lines between baseline expression
  (DepMap 24Q4) and drug AUC. Expression and AUC are both centred within lineage.
  - PRISM 19Q4 gives 1,390 drugs over a median of 430 lines.
  - For the anchor drugs, standardised PRISM, GDSC1 and GDSC2 programs are
    averaged: cisplatin, oxaliplatin, 5-FU, mitomycin C, doxorubicin, paclitaxel,
    temozolomide, cyclophosphamide.
- **Score.** `score(N, d) = −corr(s_N, r_d)` over landmark genes, where `r_d` is the
  drug-specific program. The mean program over all 1,390 drugs is projected out,
  removing the "slow-growing lines resist everything" axis. A positive score means
  the NPI moves expression toward cells that are sensitive to the drug.
- **NPI signatures.**
  - The 33 modifier signatures from `Cancer Resource/data/modifier_profiles`
    (Physiological-Fitness-Landscape, 7f72fbb).
  - Three new amino-acid deprivation signatures from GSE62673: MCF7 −cystine,
    MCF7 −methionine and PC3 0 vs 20 µM methionine, all at 24 h.

## Drug programs reproduce known resistance biology

| gene | drug | program (resistance correlation) | percentile among 1,390 drugs |
|---|---|---|---|
| ABCB1 (efflux pump) | paclitaxel / doxorubicin | +0.34 / +0.19 | 99th / 97th |
| SLC7A11 (cystine importer) | erastin | +0.21 | 100th |
| MGMT | temozolomide | +0.17 | 100th |
| SLFN11 | mitomycin C / cisplatin (PRISM) | −0.12 / −0.05 | 2nd / 16th |
| SLFN11 | cisplatin (GDSC1 / GDSC2) | −0.28 / −0.41 | |

**Reliability** (PRISM split-half r; a full-sample program is more reliable):

| drug | split-half r |
|---|---|
| doxorubicin, paclitaxel | 0.46 |
| mitomycin C | 0.34 |
| oxaliplatin | 0.28 |
| erastin | 0.27 |
| cisplatin, 5-FU | 0.14 |
| carboplatin | 0.05 |
| temozolomide, cyclophosphamide | ≈ 0 (inactive in the 5-day assay; cyclophosphamide is a prodrug) |

**PRISM vs GDSC agreement:**

| drug | r |
|---|---|
| mitomycin C | 0.59 |
| oxaliplatin (GDSC2) | 0.57 |
| doxorubicin | 0.53 |
| cisplatin | 0.10–0.26 |

## Pre-registered gates (fixed before any NPI was scored)

| gate | evidence | result |
|---|---|---|
| G1: in-vitro heat scores platinum above 5-FU and mitomycin C | tier-1: heat synergy with platinum, additive with 5-FU and MMC (Helderman 2020) | nominal pass (4/6 signatures, mean +0.011), but **not above chance**: p = 0.39 against random drug pairs, null SD 0.036 |
| G2: cystine deprivation × erastin positive and in the top 10% | tier-1: cystine withdrawal + erastin | **fail**: +0.066, 77th percentile |
| G3: methionine restriction × 5-FU and × oxaliplatin positive | tier-2 | **fail**: 3/4 positive; MCF7 × oxaliplatin is −0.135 |

Not testable: glucose restriction × metformin (metformin is in neither PRISM nor
GDSC); substrate stiffness (no signature); fasting (no cell-line signature).

Without the gates, the top-ranked pairs are dominated by obscure PRISM compounds
(KY02111, ZLN005, tiagabine). That is what noise looks like in a 36 × 1,390
screen. The full table is written to `out/cancer_npi_drug/npi_drug_scores.tsv`
(not committed) and should not be read as recommendations.

## Why it fails, and what that implies

- **Signature to sensitivity.** Baseline differences between cell lines are not
  the same thing as an acute stress response. Heat and starvation signatures at
  4–24 h are dominated by HSP, UPR and ISR programs. The best-understood NPI ×
  drug synergies act mainly after transcription:
  - heat raises platinum uptake and inactivates homologous recombination;
  - cystine withdrawal and erastin both deplete glutathione.

  An mRNA signature, especially one on 978 genes, sees these poorly.
- **Program reliability.** The platinum and 5-FU programs, the drugs that matter
  most for the gates, are the least reliable.
- **Now tested at scale.** The same scoring rule was run against DrugComb's
  740k measured drug-pair outcomes, where labels are abundant:
  no interaction signal in any of four synergy metrics
  ([docs/validation_drugcomb.md](validation_drugcomb.md)). The gate failures below
  are not a small-sample artefact.
- **Two families now fail the same known results.** Phase 7 (trained GBT) and
  this unsupervised predictor both fail. Signature-only approaches have not
  reproduced the measured NPI × drug results. The limiting input is measured
  NPI × drug outcomes, not model choice.

## Next steps that could change the verdict

1. Curate measured survival fractions for NPI alone, drug alone and combined
   (hyperthermia 41–43 °C, amino-acid and glucose deprivation, with platinum,
   5-FU, doxorubicin and paclitaxel) in the DepMap lines, and use Bliss excess
   as labels. Even 50–100 labels would allow a supervised check of any predictor.
2. Add mechanism features the transcriptome misses, then test them with the same
   gates:
   - drug uptake transporters;
   - HR/FA pathway dependency from DepMap CRISPR screens;
   - GSH and cystine dependency.
3. Keep the gates. Any future predictor must pass G1–G3 before its grid is used.

## Reproduce

```bash
# DepMap 24Q4 (figshare 27993248): Model.csv, OmicsExpressionProteinCodingGenesTPMLogp1.csv
# PRISM 19Q4 (figshare 9393293): secondary-screen-dose-response-curve-parameters.csv
# GDSC release 8.2 (ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-8.2/)
#   -> data/raw/depmap/
npi-pharma fetch-geo GSE62673 --raw-dir data/raw
python scripts/cancer_npi_drug.py --cancer-resource "<Physiological-Fitness-Landscape>/Cancer Resource"   # ~4 min first run
```
