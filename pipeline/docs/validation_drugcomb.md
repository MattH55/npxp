# Does signature composition predict measured synergy? Test on DrugComb

Run 2026-09-28 on DrugComb v1.5 (`Cancer Resource/data/drugcomb/summary_v_1_5.csv`,
1,432,351 rows; 739,964 real drug pairs over 288 cell lines and 4,268 drugs).

**Result: no. Across four synergy metrics and three increasingly strict designs,
the signature-composition scores carry no measurable information about drug
interaction.** The one apparent signal is a monotherapy-potency confound, and
even the monotherapy claim is at best marginal.

This matters because it is the check the project could not previously run. There
are almost no measured NPI × drug outcomes, but there are ~740k measured
drug × drug outcomes. The scoring rule is the same in both cases, so drug × drug
is a fair and far better powered test of the rule itself.

## What was tested

`npi_pharma.interact.score`'s reversal and complementarity, and
`npi_pharma.cancer`'s resistance-reversal, with LINCS L1000 Phase II signatures
on 949–978 landmark genes, DepMap 24Q4 expression, and PRISM/GDSC AUC.

Every test reports a **residual** correlation: the synergy minus each drug's mean
and each cell line's mean. Without that, a broadly cytotoxic drug (many synergies,
strong signature) or a fragile cell line inflates the raw correlation. The raw
numbers are given too.

| test | design | n |
|---|---|---|
| 1. pair level | `score(A,B) = mean(−corr(s_A, r_B), −corr(s_B, r_A))` vs the pair's median synergy | 2,734 pairs, 164 drugs, 159,471 measured rows |
| 2. cell-wise | npi_pharma complementarity with the cell line in the patient slot; consensus signatures | 50,069 (pair, cell line), 2,290 pairs, 165 lines |
| 3. cell-matched | same, but both signatures measured **in that same cell line** | 4,127 observations, 1,955 pairs, 9 lines |
| 4. monotherapy | reversal vs the drug's own AUC in that line | 189,924 (PRISM), 61,285 (GDSC1), 31,100 (GDSC2) |

## 1. Pair level

| synergy metric | reversal ρ raw | reversal ρ residual | p | AUC synergy vs antagonism | \|cos\| baseline ρ residual |
|---|---|---|---|---|---|
| Bliss | +0.078 | +0.017 | 0.39 | 0.617 (n+=70, n−=127) | −0.017 |
| Loewe | −0.012 | −0.029 | 0.12 | 0.455 | +0.034 |
| ZIP | +0.076 | +0.036 | 0.07 | 0.689 (n+=61, n−=5) | −0.030 |
| HSA | +0.036 | +0.005 | 0.78 | 0.491 | −0.074 |

The raw correlation of about +0.08 is almost entirely drug main effects: it falls
to +0.02 once those are removed. No metric reaches significance, and the residual
signs disagree. The `|cos|` baseline (the monotherapy-correlation feature of
eLife 2020;9:e52707) fails here too.

## 2. Cell-wise (which lines synergise?)

The question the tool is actually asked. Residual Spearman, after removing pair
and cell-line means:

| metric | complementarity | complementarity_gain | mono_reversal | reverse_combo | p (complementarity) |
|---|---|---|---|---|---|
| Bliss | −0.007 | −0.005 | +0.005 | +0.005 | 0.20 |
| Loewe | −0.001 | −0.000 | −0.004 | −0.004 | 1.00 |
| ZIP | −0.003 | −0.003 | −0.008 | −0.008 | 0.53 |
| HSA | −0.016 | −0.012 | +0.037 | +0.036 | 0.10 |

Zero, on 50,069 measured observations.

## 3. Cell-matched (the strictest design)

Both drugs profiled by LINCS in the same cell line the synergy was measured in
(HT29, A375, MCF7, PC3, A549, HS578T, MDA-MB-231, LNCaP, BT20):

| metric | complementarity | mono_reversal | reverse_combo |
|---|---|---|---|
| Bliss | −0.004 (p = 0.80) | **−0.218** | **−0.221** |
| Loewe | +0.008 (p = 0.62) | −0.144 | −0.147 |
| ZIP | −0.026 (p = 0.10) | **−0.260** | **−0.263** |
| HSA | +0.075 (p = 0.0005) | **+0.246** | **+0.251** |

Complementarity, the interaction term, is still zero. But `mono_reversal` and
`reverse_combo` are large — **and flip sign between metrics**. That is diagnostic,
not encouraging: when both drugs are potent in a line, Bliss and ZIP excess fall
toward their ceiling while HSA excess rises. A feature tracking "both drugs work
here" produces exactly this pattern. So the score is reading monotherapy potency,
not interaction. The one significant complementarity result (HSA, +0.075) is the
same confound leaking through the metric most coupled to potency.

## 4. So does it predict monotherapy potency?

Barely, and not consistently:

| source | n | ρ raw | ρ residual | p | per-drug median | fraction in the predicted direction |
|---|---|---|---|---|---|---|
| PRISM | 189,924 | −0.025 | **+0.005** | 0.027 | +0.004 | 47% |
| GDSC1 | 61,285 | −0.089 | −0.022 | 0.003 | −0.031 | 60% |
| GDSC2 | 31,100 | −0.042 | −0.021 | 0.003 | −0.015 | 64% |

Negative is the predicted direction (more reversal, lower AUC, more sensitive).
GDSC agrees with the prediction and PRISM contradicts it. At n = 190k, even
ρ = 0.005 is "significant"; the effect size, not the p-value, is what matters
here, and |ρ| ≤ 0.022 is not usable.

## What this settles

- **Three model families have now failed on measured data.** Phase 7's trained
  pathway-feature GBT (see `analysis/phase7_audit/`), the unsupervised
  resistance-reversal predictor (`docs/cancer_npi_drug_predictor.md`), and the
  npi_pharma complementarity score all fail to reproduce measured interaction.
  The common factor is the input, not the architecture: an L1000 signature of an
  agent plus a baseline expression profile does not determine how two agents
  interact.
- **This is consistent with the mechanisms.** Heat potentiates platinum by raising
  drug uptake and disabling homologous recombination; cystine withdrawal and
  erastin converge on glutathione depletion. These act on transport, protein
  state and metabolite pools, largely downstream of or parallel to mRNA.
- **The npi_pharma `complementarity` score should not be read as a synergy
  prediction**, including in the NPI × drug reports. It remains what its
  disclaimer says, a mechanism-only plausibility ranking, and now with a measured
  bound: no detectable relation to synergy in 50k drug pair observations.

## What would be worth trying next

1. **Features beyond the transcriptome**, tested against these same DrugComb
   labels before any NPI use: drug influx/efflux transporter expression and
   copy number, DepMap CRISPR dependency on the pathway the partner drug damages
   (HR, FA, base excision repair), and GSH/cystine metabolism state. The tests
   in this directory make that a one-command check.
2. **Measured NPI × drug outcomes.** Even 50–100 curated triplets (NPI alone,
   drug alone, combined) in DepMap cell lines would allow direct supervision
   instead of transfer from drug pairs.
3. **Keep the gates.** `scripts/cancer_npi_drug.py` holds the pre-registered
   tier-1 gates; these four tests are the drug × drug counterpart. Any new
   predictor should clear both before its ranking is reported.

## Reproduce

```bash
# DrugComb v1.5 is tracked with Git LFS in MattH55/Physiological-Fitness-Landscape
# (Cancer Resource/data/drugcomb/summary_v_1_5.csv, 1.42 GB): git lfs pull --include=...
DC="<checkout>/Cancer Resource/data/drugcomb/summary_v_1_5.csv"
npi-pharma ingest-lincs --gctx ... --gene-space landmark --perts "<matched drugs>" --out drugs.parquet
python scripts/validate_reversal_drugcomb.py  --drugcomb "$DC" --lincs drugs.parquet
python scripts/validate_cellwise_drugcomb.py  --drugcomb "$DC" --lincs drugs.parquet
python scripts/validate_cellmatched_drugcomb.py --drugcomb "$DC"
python scripts/validate_monotherapy.py --lincs drugs.parquet
```
