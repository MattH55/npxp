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

## Why this matters for the project

Three model families have now failed to predict interaction from signatures. The
diagnosis throughout was that the missing ingredient is measured interaction data,
not a better model. This is a route to that data:

1. Scale the scan — 16% of screened series were factorial, and only 70 were screened.
2. Fix the expression-mapping failures, which cost 6 of 11 series.
3. With enough of them, ask directly whether any single-agent feature predicts the
   measured interaction contrast — the experiment that was impossible before.

## Reproduce

```bash
python scripts/find_combination_series.py --max-candidates 60 --max-fetch 70
python scripts/build_interaction_signatures.py
```
