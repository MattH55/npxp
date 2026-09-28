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
| 5. positive control | measured same-experiment monotherapy response vs synergy | 739,964 rows, 73,167 pairs, 288 lines |

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

## 5. Positive control: is the residual bar simply unreachable?

No. The identical residual test, run on the feature the literature calls the
strongest single predictor of synergy (each agent's own response, measured in
the same DrugComb experiment):

| metric | ri_min | ri_mean | ri_max | css_mean |
|---|---|---|---|---|
| Bliss | −0.257 | **−0.302** | −0.257 | +0.083 |
| Loewe | −0.016 | −0.137 | −0.188 | +0.059 |
| ZIP | −0.223 | −0.243 | −0.203 | +0.076 |
| HSA | −0.153 | −0.224 | −0.221 | +0.056 |

Residual Spearman on 739,964 rows, 73,167 pairs, 288 cell lines. So the test
detects real structure at |ρ| ≈ 0.2–0.3 where it exists, an order of magnitude
above the ≤ 0.02 the signature scores reach. The nulls above are about the
features, not the bar.

Two honest notes on this control. The sign says more single-agent inhibition
gives *less* Bliss and ZIP excess, which is largely the metric's ceiling: if
each drug alone already kills the cells, there is little headroom for excess.
So it is a strong statistical relation rather than proof that monotherapy
response predicts biological synergy. That is consistent with the DREAM
challenge, where a monotherapy-only predictor matched the average submitted
model (Menden et al., Nat Commun 2019, doi:10.1038/s41467-019-09799-2).

## How this sits in the literature

Drug-synergy prediction from cell-line features is a mature field, and these
results are consistent with it rather than contradicting it.

- The DREAM challenge (160 teams, 11,576 experiments, 85 cell lines) found
  synergy predictable to replicate-level accuracy for >60% of combinations, that
  winning methods needed prior drug-target knowledge, and that 20% of
  combinations were poorly predicted by every method.
- A 2024 review of newer methods (arXiv:2404.02484) reports that the best models
  solve scenarios involving *known* drugs and cell lines (AUROC up to 0.98,
  Pearson up to 0.89), while "scenarios involving new drugs or cell lines still
  fall short of an accurate prediction level", with leave-drug-out significantly
  harder than leave-pair-out.
- Most relevant: a 2025 study (Front Pharmacol, PMC12310685) builds "Drug
  Resistance Signature" features in nearly the way `npi_pharma.cancer` does, from
  LINCS signatures plus GDSC-defined resistant and sensitive lines, and reports
  large gains (Pearson 0.72 vs 0.68; MSE 92 vs 346). But it evaluates with
  stratified random 5-fold cross-validation, with no leave-drug-out or
  leave-cell-line-out split and no main-effect ablation. Under a random split the
  same drugs, pairs and cell lines appear in training and test, so drug and
  cell-line main effects, which this project's audits show dominate, are
  learnable rather than held out.

The contribution here is therefore not a new finding about drug pairs. It is the
effect size after main effects are removed, on the specific scoring rule this
project uses, in the regime this project needs: a novel agent in a given cell
line. That regime is the one the field itself reports as unsolved.

## What this settles

- **The bar is fair, and the features fail it.** Monotherapy response clears the
  same test at |ρ| ≈ 0.3; the signature scores reach ≤ 0.02.
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
python scripts/validate_positive_control.py --drugcomb "$DC"
```
