# Applying the literature's finding: rank by efficacy, from single-agent response

Run 2026-09-28. Code: `src/npi_pharma/cancer/efficacy.py`,
`scripts/efficacy_from_monotherapy.py`, `scripts/rank_npi_drug_efficacy.py`,
`configs/npi_monotherapy.yaml`.

The synergy-prediction literature converges on two points: each agent's own
response in the cell line is the strongest single predictor of a combination's
outcome (DREAM challenge, Menden et al., Nat Commun 2019), and predicting
excess-over-additivity for unseen agents remains unsolved (review,
arXiv:2404.02484). This project's own tests agree — signature composition reaches
|ρ| ≤ 0.02 on measured synergy (`validation_drugcomb.md`).

So stop predicting synergy and predict the quantity that is predictable.

## Efficacy and synergy are nearly independent

In DrugComb, across 738,619 drug-pair experiments:

| | correlation with combination efficacy (CSS) |
|---|---|
| agents' own mean relative inhibition | **+0.68** |
| Bliss synergy | +0.06 |

Combination efficacy is how much the combination actually inhibits. Synergy is
how much it exceeds additivity. They measure different things, and it is efficacy
that decides whether a combination kills the cells.

## How well each is predictable in an unseen cell line

Leave-cell-lines-out (GroupKFold by cell line; every test cell line unseen),
738,619 rows, 73,114 pairs, 236 cell lines:

| target | predictor | Spearman | Pearson | RMSE | within-cell-line median ρ |
|---|---|---|---|---|---|
| **efficacy (CSS)**, sd 25.5 | pair mean elsewhere | 0.708 | 0.682 | 18.68 | 0.619 |
| | single-agent response | 0.706 | 0.724 | 17.58 | 0.725 |
| | both | **0.809** | **0.815** | **14.77** | **0.776** |
| **Bliss synergy**, sd 47.4 | pair mean elsewhere | 0.226 | 0.041 | 47.95 | 0.326 |
| | single-agent response | 0.458 | 0.299 | 45.19 | 0.335 |
| | both | 0.535 | 0.289 | 45.41 | 0.443 |

Efficacy: RMSE 14.77 against a target sd of 25.5, so about two thirds of the
variance (R² ≈ 0.66). Synergy: RMSE 45.4 against sd 47.4, about 9%.

## The parameter-free form, which is what transfers to NPIs

A trained model needs training examples of the agents. NPIs have none. Bliss
independence of the two measured single-agent effects, `1 − (1−a)(1−b)`, needs no
fitting at all, and on the same 739,964 rows:

| predictor | Spearman with measured efficacy | within-cell-line median ρ |
|---|---|---|
| Bliss expectation, `1 − (1−a)(1−b)` | 0.696 | **0.752** |
| sum of single-agent inhibitions | 0.690 | |
| max of the two | 0.688 | |
| trained GBT (for reference) | 0.809 | 0.776 |

The parameter-free form gives up about 0.11 Spearman and buys the ability to
score an agent with no training data. That is the trade this project needs.

## Applied to NPI × drug

`scripts/rank_npi_drug_efficacy.py` ranks NPI × drug per cell line with the NPI's
curated single-agent inhibition and the drug's PRISM/GDSC AUC. The blocker is
now explicit and is data, not method: `configs/npi_monotherapy.yaml` holds
**3 curated entries covering 2 cell lines**, all hyperthermia at 43 °C/60 min
from Helderman et al. 2020:

| NPI | cell line | assay | inhibition |
|---|---|---|---|
| 43 °C / 60 min | RKO | clonogenic, 10 d | 0.81 ("only 19% of the cells surviving") |
| 43 °C / 60 min | RKO | short-term viability | 0.30 ("a 30% reduction") |
| 43 °C / 60 min | HCT116 | short-term viability | 0.0 ("does not show any differences") |

With the drug side from PRISM and GDSC this ranks 3,374 NPI × drug rows over
1,680 drugs in those 2 lines. Predicted combined inhibition, HIPEC drugs:

| drug | drug alone, RKO | RKO + heat | drug alone, HCT116 | HCT116 + heat | measured TER at 43 °C |
|---|---|---|---|---|---|
| mitomycin C | 0.80 | 0.96 | 0.49 | 0.49 | not increased |
| 5-fluorouracil | 0.41 | 0.89 | 0.45 | 0.45 | not increased |
| oxaliplatin | 0.22 | 0.85 | 0.17 | 0.17 | 3.3 (RKO), 3.1 (HCT116) |
| cisplatin | 0.13 | 0.84 | 0.09 | 0.09 | 3.5 (RKO), 2.8 (HCT116) |
| carboplatin | 0.00 | 0.81 | 0.06 | 0.06 | 7.2 (RKO) |

### Read this carefully — three honest limits

1. **It ranks total kill, not synergy, and the two disagree here.** The top
   predicted combination in RKO is heat + mitomycin C (0.96), yet mitomycin C is
   one of the two drugs Helderman reports as having *no* thermal synergy. There is
   no contradiction: mitomycin C is simply potent on its own (0.80), so the
   combination kills a lot without any interaction. If what you want is
   interaction, this ranking does not supply it — and nothing in this repository
   currently does.
2. **Total kill may be the wrong target for HIPEC.** Heat kills normal tissue too.
   The reason to prefer a thermally synergistic drug is tumour selectivity, which a
   total-kill ranking ignores. Helderman ran normal colon organoids for exactly
   this reason. A selectivity target needs paired normal-tissue measurements.
3. **Two cell lines cannot validate anything.** The ranking puts every combination
   higher in RKO than HCT116, because heat alone contributes 0.81 versus 0.0. The
   measured TERs are also higher in RKO (cisplatin 3.5 vs 2.8), which is
   directionally consistent — but with n = 2 that is an observation, not evidence.

## Curation attempted, and why it is the gate

Five systematic literature passes (PMC open-access full text plus PubMed
abstracts) added **nothing** to the 3 entries below, and a DepMap CRISPR
dependency proxy was tested and rejected. See
[npi_monotherapy_search_log.md](npi_monotherapy_search_log.md) for what was
searched, why keyword search cannot reach this data, and the calibration finding
that a per-cell-line functional measurement tops out near +0.2 to +0.47 Spearman
for single-agent sensitivity even in the most favourable case.

## What to curate next

`configs/npi_monotherapy.yaml` lists the wanted measurements in priority order.
Each needs a published surviving fraction or viability loss for the NPI **alone**
in a named cell line; the thermal enhancement ratios already curated are a
different quantity and cannot substitute.

1. Hyperthermia 41/42/43 °C / 60 min in the other 11 curated cell lines. The
   classic thermal-dose literature reports these as survival vs CEM43.
2. Glucose restriction (5% medium glucose, 96 h) in T47D, MCF7, DU145.
3. Cystine and methionine withdrawal in MDA-MB-231 and HCT116.
4. Serum starvation in LoVo; acidosis and hypoxia in U-87 MG.

With those, the same command ranks the full NPI × drug × cell-line grid on
measured inputs rather than transcriptomic inference.

## Reproduce

```bash
DC="<checkout>/Cancer Resource/data/drugcomb/summary_v_1_5.csv"
python scripts/efficacy_from_monotherapy.py --drugcomb "$DC"     # ~90 s
python scripts/rank_npi_drug_efficacy.py                         # needs data/raw/depmap
```
