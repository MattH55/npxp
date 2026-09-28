# Phase 7 audit: can the predictions rank NPI × drug combinations?

Audit of `Cancer Resource/data/lincs/phase7_predictions.json` in
MattH55/Physiological-Fitness-Landscape (commit 7f72fbb): 3,861 predictions,
33 modifier signatures × 9 drugs × 13 cell lines, from the pathway-feature GBT.
Reproduce with:

```bash
python analysis/phase7_audit/audit_phase7.py --cancer-resource "<checkout>/Cancer Resource" --out results.json
```

**Verdict:** Phase 7's RMSE win over Bliss (0.108–0.109 vs 0.111) does not carry
over to ranking combinations. Its predictions fail every check below that it was
not trained on, including the one direct comparison with known biology. Keep them
`validated: false`; do not prioritise combinations from them.

## Checks

| check | expected if the ranking carries signal | result |
|---|---|---|
| Phase 3 negative controls (weakest mimetic anchoring) vs the other 30 modifiers | AUC of a real modifier outscoring a control > 0.5 | **0.41**; fasting 40 h (a control) ranks 3rd of 33 |
| predictions above an approximate 90% half-width (1.645 × drug-out RMSE = 0.178) | a meaningful set | **8 of 3,861**, and drug-out coverage is already under-nominal (0.823) |
| variance share | mostly interaction | modifier 0.28, cell line 0.17, drug 0.13, interaction residual 0.42 |
| GSE85620 arms: same 10-week strength programme, placebo drink vs cold-water immersion | similar scores | ranked **2nd** (mean 0.058) vs **16th** (0.015) |
| curated heat × cisplatin enhancement ratios (below) | heat scores positive and above other modifiers | **reversed** for heat applied to cancer cells |

## Directional anchor check

Three curated tier-1 values (`tests/test_seed.py` `REAL_COMBINED_EFFECT_METRICS`)
have a drug and cell line inside the Phase 7 grid. All three are synergy. The
metric differs from the model's target (enhancement ratio vs Bliss excess of
survival fraction), so only direction and rank are compared. The carboplatin,
oxaliplatin and stiffness anchors fall outside the grid.

| cell line (cisplatin) | literature enhancement ratio | heat applied to cancer cells in vitro (6 signatures): mean pred, positive, mean rank / 33 | heat in muscle / blood in vivo (3): mean pred, positive, mean rank / 33 | all 33 modifiers, mean |
|---|---|---|---|---|
| HeLa | 3.10 (Kusumoto 1993, 42.8 °C) | −0.076, 0/6, 28.7 | +0.039, 3/3, 8.3 | −0.005 |
| RKO | 3.5 (Helderman 2020, 43 °C) | −0.001, 2/6, 25.2 | +0.091, 3/3, 11.7 | +0.050 |
| HCT116 | 2.8 (Helderman 2020, 43 °C) | +0.005, 2/6, 26.7 | +0.102, 3/3, 11.7 | +0.064 |

The in-vitro heat signatures (GSE48398 45 °C/30 min in MCF7, MDA-MB-231,
MDA-MB-468 and MCF10A; GSE10043 U937 41 °C; GSE75127 HSC-3 44 °C) are the
closest analogue of the anchor experiments, yet they score near the bottom.
The in-vivo muscle and blood heat signatures score positive, but they share
little with heating a tumour. This matches the domain gap in the build order:
the model is trained on drug × drug pairs, and most modifier signatures are
from non-tumour tissue.

## What the model ranks highest (for reference, not as recommendations)

Mean predicted Bliss excess over the 13 cell lines:

| modifier | drug | mean | range |
|---|---|---|---|
| muscle heat (GSE82323) | paclitaxel | 0.104 | 0.049–0.183 |
| strength training, placebo drink (GSE85620) | paclitaxel | 0.103 | 0.021–0.192 |
| cold exposure 10 d (GSE156248) | erastin | 0.100 | 0.061–0.147 |
| exhaustive exercise, WBC (GSE3606) | cisplatin | 0.097 | 0.072–0.123 |
| cold exposure 10 d | 5-fluorouracil | 0.091 | 0.063–0.118 |
| strength training, placebo drink | cisplatin | 0.091 | 0.041–0.189 |

## Next

- Make the heat × platinum direction a gate. Any future model must score in-vitro
  heat with cisplatin positive in HeLa, RKO and HCT116 before its grid is used.
- Extend the anchors inside the grid: curate in-vitro hyperthermia with doxorubicin,
  paclitaxel and mitomycin C survival fractions (drug alone, heat alone, combined)
  for the 13 cell lines. That gives Bliss excess directly, in the model's own metric.
- Put tumour- or cell-line-derived modifier signatures ahead of muscle and blood ones.
