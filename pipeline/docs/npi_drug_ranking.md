# NPI × drug ranking across every drug signature, with reliability flags

Run 2026-09-29. Code: `scripts/rank_npi_drug_all.py`.
Output: `out/npi_drug_all/npi_drug_all.tsv` (24,906 rows).

14 human in-vitro cancer-cell-line NPI signatures scored against **1,779 drug
signatures**. Weak signatures are included and labelled rather than withheld, so
every row can be read with its confidence attached.

## What the score means

`cos > 0`: the NPI moves expression the way the drug does — a **mimic**, so pairing
them risks duplicating one mechanism rather than adding two. `cos < 0`: it opposes.

**Similarity is not synergy.** This project measured that neither similarity nor
orthogonality predicts it ([validation_drugcomb.md](validation_drugcomb.md)). Use
these as mechanism hypotheses, never as combination recommendations.

## Reliability bands

| band | n | basis |
|---|---|---|
| high | 1,763 | LINCS consensus: a median over ~7 cell lines per compound |
| medium | 2 | GEO consensus, cross-series agreement ≥ 0.30, from independent studies |
| low | 3 | GEO consensus, agreement 0.10–0.30 |
| very low | 11 | GEO consensus below 0.10, or a single-series signature |

The bands come from measurement, not judgement. Validating the GEO pipeline against
LINCS on 10 drugs gave 0.39 median agreement, and cross-series agreement predicts
that at r = 0.70 ([drug_consensus.md](drug_consensus.md)). Single-cell-line
signatures are capped at "very low" because the same drug in two cell lines agrees
at only 0.18, while two *different* drugs in one cell line agree at 0.71
([npi_drug_retrieval.md](npi_drug_retrieval.md)).

Two independence rules cap a consensus regardless of its number: arms from one
accession are not independent evidence, and **consecutive accessions are almost
always companion series of one submission** — GSE59296 and GSE59297 are a single
topotecan study split by platform, which is why topotecan is "very low" despite
having two "series".

| GEO consensus | agreement | band |
|---|---|---|
| erastin | 0.683 | medium |
| actinomycin D | 0.347 | medium |
| 5-fluorouracil | 0.219 | low |
| cyclophosphamide | 0.120 | low |
| melphalan | 0.107 | low |
| topotecan, cisplatin | 0.016 | very low |
| 5-azacytidine | 0.008 | very low |
| carboplatin | −0.040 | very low |
| methotrexate | −0.049 | very low |

## Top high-reliability match per NPI

| NPI | cell line | closest drug (LINCS consensus) | cosine |
|---|---|---|---|
| serum starvation 96 h | LoVo | IKK-2 inhibitor V | 0.452 |
| glucose deprivation | T47D | palbociclib | 0.522 |
| hyperthermia 44 °C | HSC-3 | SMER-3 | 0.374 |
| heat 45 °C | MDA-MB-231 | camptothecin | 0.349 |
| heat 45 °C | MCF10A | ruxolitinib | 0.330 |
| glucose deprivation | MCF7 | cepharanthine | 0.327 |
| heat 45 °C | MCF7 | SMER-3 | 0.285 |
| mild hyperthermia 41 °C | U937 | radicicol | 0.254 |
| β-hydroxybutyrate 10 mM | MCF7 | propylthiouracil | 0.243 |
| heat 45 °C | MDA-MB-468 | eplerenone | 0.194 |
| β-hydroxybutyrate 25 mM | T47D | leflunomide | 0.165 |
| cystine deprivation | MCF7 | RAF-265 | 0.057 |
| methionine deprivation | MCF7 | QL-XII-47 | 0.046 |
| methionine deprivation | PC3 | VX-745 | 0.025 |

Two of these are mechanistically coherent rather than arbitrary: **heat in
MDA-MB-231 resembling camptothecin** (a topoisomerase poison — heat also produces
DNA damage and replication stress), and **heat in U937 resembling radicicol** (an
HSP90 inhibitor, which triggers the heat-shock response directly). Glucose
deprivation in T47D resembling palbociclib (CDK4/6 inhibition) fits both causing
G1 arrest.

The three amino-acid-deprivation NPIs top out at 0.03–0.06 — essentially nothing
resembles them in a 1,763-compound panel. That is a real negative result, not a
failure to look.

## Where the curated drugs land

Best row per curated drug, across all NPIs:

| drug | best NPI | cosine | reliability |
|---|---|---|---|
| erastin | glucose deprivation, T47D | 0.412 | very low (single series) |
| 5-fluorouracil | serum starvation, LoVo | 0.356 | very low (single series) |
| carboplatin | heat, MCF10A | 0.334 | very low (single series) |
| **doxorubicin** | serum starvation, LoVo | **0.285** | **high** |
| oxaliplatin | serum starvation, LoVo | 0.245 | very low (single series) |
| **mitomycin C** | serum starvation, LoVo | **0.185** | **high** |
| cisplatin | methionine deprivation, MCF7 | 0.173 | very low |
| **paclitaxel** | heat, MCF7 | **0.168** | **high** |
| **temozolomide** | heat, MCF10A | **0.100** | **high** |
| **metformin** | heat, MCF7 | **0.075** | **high** |

Only five curated drugs have a high-reliability signature, and the platinums are
not among them — the reason the platinum consensus work matters. Read the "very
low" rows as provisional.

## Honest reading

- The LoVo serum-starvation signature tops several lists. Serum starvation is a
  broad stress, so its resemblance to many drugs may reflect a generic stress
  program rather than specific mechanism. Its high scores should be discounted.
- Every similarity here is cross-cell-line except where noted, and cell line
  dominates drug identity in single-series signatures. The LINCS consensus rows
  are the ones least affected.
- Nothing here licenses a combination recommendation.

## Reproduce

```bash
python scripts/rank_npi_drug_all.py --cancer-resource "<...>/Cancer Resource"
```
