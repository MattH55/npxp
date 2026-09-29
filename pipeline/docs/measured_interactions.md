# Measured interaction signatures from factorial GEO designs

Run 2026-09-29. Code: `src/npi_pharma/interact/measured.py`,
`scripts/find_combination_series.py`, `scripts/build_interaction_signatures.py`.

Every interaction score elsewhere in this package is *inferred* from single-agent
signatures, and that was measured not to work: signature composition carries no
detectable information about synergy across 50,000 drug-pair observations
([validation_drugcomb.md](validation_drugcomb.md)). A GEO series that profiles
control, agent A, agent B **and** the A+B combination is different in kind,
because the interaction can be computed instead of predicted:

    I(A, B) = (combo − control) − [(A − control) + (B − control)]
            = combo − A − B + control

on a log-expression scale. That is the standard two-factor interaction contrast and
the transcriptional analogue of a Bliss excess: genes with large positive I are
induced by the combination beyond what the two agents do separately, and large
negative I means suppressed beyond additivity.

## These series exist, and in useful numbers

GEO holds roughly 988 human expression series mentioning synergy and 602 naming a
drug combination explicitly. Screening 70 of them found **11 with a true factorial
design** — a 16% hit rate — so the full corpus plausibly holds well over a hundred.

| series | agents | combination arm |
|---|---|---|
| GSE338229 | **cisplatin** + nelfinavir | yes |
| GSE311210 | palbociclib + fulvestrant | yes |
| GSE336734, GSE336729 | abemaciclib + fulvestrant | yes |
| GSE270318 | venetoclax + ERK inhibitor | yes |
| GSE325471 | dalpiciclib + enzalutamide | yes |
| GSE262443 | A51 + ceritinib | yes |
| GSE304294 | IMMU132 + IACS010759 | yes |
| GSE315147 | ACT001 + Stattic | yes |
| GSE294096 | VTP50469 + WM119 | yes |
| GSE330930 | B7H3 + copper | yes |

Several are real clinical pairings (palbociclib or abemaciclib with fulvestrant is
standard in HR-positive breast cancer), and one uses cisplatin, a drug carrying much
of this project's tier-1 evidence.

## Five built so far

| series | agents | genes | sd of interaction | arms |
|---|---|---|---|---|
| GSE338229 | cisplatin + nelfinavir | 14,184 | 0.322 | 3/3/3/3 |
| GSE270318 | venetoclax + ERK inhibitor | 12,826 | 0.442 | ok |
| GSE262443 | A51 + ceritinib | 13,542 | 0.589 | ok |
| GSE304294 | IMMU132 + IACS010759 | 16,918 | 0.545 | small_n |
| GSE336734 | abemaciclib + fulvestrant | 13,084 | 0.654 | small_n |

The remaining six failed on expression availability, not design: their supplementary
tables could not be mapped to samples through GEO's labels, and the pipeline refuses
to guess a column.

For GSE338229 (cisplatin + nelfinavir in the cisplatin-resistant ovarian line
IGROV-CP), the genes furthest below additivity include CEACAM1, COL17A1, DGKA and
MUC2 — adhesion and membrane genes suppressed by the combination beyond what either
drug does alone.

## What this is, and is not

- **Observed, not predicted.** This is the only interaction data in the repository
  that is measured rather than assumed.
- **Noisier than a main effect.** Four group means enter the contrast, so its
  standard error combines all four arms. Two of the five have fewer than three
  replicates per arm and are flagged `small_n`.
- **One context each.** Each signature is one cell line, one dose pair, one
  timepoint. It does not generalise on its own, exactly as the cell-line dominance
  result would predict ([npi_drug_retrieval.md](npi_drug_retrieval.md)).
- **Not an NPI result.** These are drug × drug. Their value here is as ground truth
  a predictor could finally be tested against, and as a template: an NPI × drug
  factorial series (heat or nutrient restriction with a drug, all four arms) would
  give the same quantity for the question this project actually asks.

## Item 2 is done: 10 of 11 series now build

The expression-mapping failures cost 6 of 11 series. Five of the six are recovered,
and the corpus is now **10 measured interaction signatures**. The sixth, GSE330930,
publishes only a 10x `feature_reference.csv` — there is no expression matrix to map,
and it is correctly refused rather than forced.

The failures were all the same thing: RNA-seq count tables whose columns are *arm
labels* (`CAMA-1_DMSO_1`, `Veh_R`, `DMSO_R1`, `SNB19ACT`) rather than GSM names, so
no single label matched a column outright. `match_columns_by_assignment` in
`scripts/build_drug_consensus.py` solves all the samples at once instead of one at a
time, which is what makes those tables tractable: a column has to beat every *other
sample* as well as every other column, so a control arm can be identified purely by
being the column left over (`SNB19` against `SNB19ACT`/`SNB19Sta`/`SNB19Com`).

Three conditions must hold or the table is refused, and the third is the one that
keeps this from being guesswork: forbidding any assigned pair must make the total
strictly worse. A table whose labels are genuinely ambiguous produces a tie, and a
tie is a refusal. A single-sample request is refused outright, because with one
sample the joint constraint that justifies the method does not exist.

One further bug surfaced: GSE311210 indexes genes by versioned Ensembl IDs, and left
untranslated it shared *zero* genes with every other series — so it dropped silently
out of every comparison instead of failing loudly. `load_ensembl_map` now translates
them.

## Item 3, run early: the interaction contrast does not reproduce

Two series in the corpus profile the same combination independently — abemaciclib +
fulvestrant, in MCF7 (GSE336734) and CAMA1 (GSE336729). That allows the question to
be asked directly, and it comes with its own positive control, because the same four
arms of the same two series also give each agent's main effect. Every step of the
pipeline is shared, so a split between them is specific to the contrast:

| contrast | cross-series cosine |
|---|---|
| combination main effect (`combo − control`) | **+0.632** |
| abemaciclib main effect | **+0.512** |
| fulvestrant main effect | **+0.194** |
| **interaction** (`combo − a − b + control`) | **−0.057** |

Over 11,403 shared genes. The main effects reproduce — the combination's at 0.63,
well above the 0.18 that the same drug in two different cell lines manages in this
project's own measurements (`npi_drug_retrieval.md`). The interaction does not: at
−0.06 it is indistinguishable from noise, and *below* the median unrelated pair in
this corpus (+0.037). Meanwhile two mechanistically unrelated combinations reach
+0.347, which is the same pathology seen throughout this project — the cell line
speaks louder than the perturbation.

This is what the arithmetic predicts. The interaction is a difference of four group
means, so its variance is the sum of all four; at n=2 and n=3 per arm it is the
noisiest quantity the design can produce, and the main effects are the least noisy.
Being *measured* rather than inferred does not rescue it.

**Flagging, not withholding.** The 10 signatures are built and stored, with their
`quality_flag` and per-arm counts attached. What this result forbids is treating any
one of them as an estimate of that combination's interaction. What it does not
forbid is using them where the noise is averaged over — which needs many more
series, not better contrasts.

### What this changes

It reorders the remaining work. Item 3 — asking whether a single-agent feature
predicts the measured interaction — cannot be run against a target that does not
reproduce, because the ceiling on any such model is the target's own reliability,
and here that is zero. Scaling the scan (item 1) is now the prerequisite rather than
an improvement: the route forward is a per-combination *consensus* over independent
series, exactly as the drug panel needed (`drug_consensus.md`), and that needs
several series per combination where the corpus currently has two for one pairing.

The honest caveat: this is **one** independent pair, at n=2 and n=3, in two different
cell lines. It establishes that the contrast can fail this badly while its own main
effects reproduce; it cannot say how often. A second replicate pair would settle
whether −0.06 is typical or unlucky, and `REPLICATE_PAIRS` in the validation script
is where one gets added.

## The main effects are worth keeping, and they are free

The interaction from these series does not reproduce, but the main effects from the
*same four arms* do (0.512, 0.632 above). `scripts/drug_signatures_from_factorial.py`
takes them: each series gives `a - control` and `b - control`, built by the same
constructor the drug panel uses, so they drop straight into the same comparisons.

**20 single-agent signatures from 10 series**, no new downloads:

| | |
|---|---|
| cisplatin, nelfinavir | GSE338229, IGROV-1/CP |
| abemaciclib, fulvestrant | GSE336734 (MCF7) and GSE336729 (CAMA1) |
| palbociclib, fulvestrant | GSE311210, ER+ breast PDX |
| dalpiciclib, enzalutamide | GSE325471, prostate |
| ABT199, ERKi | GSE270318, OCI-AML3 |
| A51, ceritinib | GSE262443, CLB-GA |
| ACT001, stattic | GSE315147, SNB19 |
| VTP50469, WM119 | GSE294096 |
| IMMU132, IACS010759 | GSE304294, KYSE30 |

Fulvestrant now has three independent series and abemaciclib two, which is the start
of the per-drug consensus the panel needs.

These are single-series, single-cell-line signatures, which this project measured to
be dominated by the cell line rather than the drug (0.71 between two *different*
drugs in one line, against 0.18 for the same drug in two lines). So they carry the
same **"very low"** band as every other single-series signature and are useful as
consensus inputs, not on their own.

### One of them is flagged, not dropped

GSE338229's line is **IGROV-1/CP** — IGROV-1 selected for cisplatin resistance. The
response of a resistant line to the drug it resists is not that drug's response, so
this signature answers a different question from the one the panel asks. The
consensus builder rejects resistant *arms* (`RESISTANCE_VALUE`) but says nothing
about the cell line, which is annotated elsewhere in the series.
`resistant_line_flag` catches the naming conventions for derivative lines
(`/CP`, `/CP70`, `-R`, `CisR`, `MDR`, spelled-out "resistant") and marks the
signature `ok,resistant_cell_line` rather than discarding it. Its false-positive risk
is the reason it is a flag: `IGROV-1`, `OCI-AML3`, `CLB-GA`, `SNB19` and `CAMA1` must
all come back clean, and the tests assert that.

This matters beyond one signature: it is the only cisplatin signature the factorial
corpus yields, and cisplatin is one of the platinums the whole consensus re-run
exists to reach.

## Why this matters for the project

Three model families have now failed to predict interaction from signatures. The
diagnosis throughout was that the missing ingredient is measured interaction data,
not a better model. This is a route to that data, and the finding above sharpens
what "enough of it" means:

1. Scale the scan — 16% of screened series were factorial, and only 70 were screened.
2. ~~Fix the expression-mapping failures~~ — done, 10 of 11 series now build.
3. Reach several independent series *per combination*, so a consensus interaction can
   be taken. Only then can item 3 — does any single-agent feature predict the measured
   interaction — be asked against a target reliable enough to be predicted.

## Reproduce

```bash
python scripts/find_combination_series.py --max-candidates 60 --max-fetch 70
python scripts/build_interaction_signatures.py
python scripts/validate_interaction_reproducibility.py
python scripts/drug_signatures_from_factorial.py
```
