# NPI → drug retrieval, with drug signatures broadened beyond LINCS

Run 2026-09-28. Code: `scripts/find_drug_datasets.py`,
`src/npi_pharma/ingest/drug_geo.py`, `configs/drug_signatures_geo.yaml`,
`scripts/npi_drug_retrieval.py`.

Retrieval is the one method in this project that held up against independent data
(GSE95640 LCD recovered a published analysis's drug hits at AUC 0.75, p = 0.0005).
It answers which drugs an NPI's transcriptional response *resembles* — not whether
they synergise, which signatures do not predict
([validation_drugcomb.md](validation_drugcomb.md)).

**Headline finding: a drug's signature depends more on the cell line it was
measured in than on the drug.** That limits cross-cell-line retrieval and explains
why the GSE95640 result worked with LINCS *consensus* signatures.

## LINCS alone is not enough

LINCS Phase II (GSE70138) carries only 5 of the 13 curated drugs — mitomycin C,
doxorubicin, paclitaxel, temozolomide, metformin — and **no platinum agent**,
no 5-fluorouracil, erastin, triapine, lomustine or cyclophosphamide. The platinums
carry most of this project's tier-1 evidence, so they had to come from GEO.

`scripts/find_drug_datasets.py` searched GEO DataSets for human expression series
naming each drug, scored them for being acute treated-vs-control cell-line
experiments, and checked whether a processed expression table is published (a
series whose only supplement is `_RAW.tar` needs array or fastq processing). Of 43
shortlisted series, 33 had processed data.

## Six new drug signatures, all from acute treatment in human cancer lines

| drug | series | cell line | design | n |
|---|---|---|---|---|
| cisplatin | GSE319411 | **MDA-MB-231** (curated) | 1 µM, 48 h vs control medium | 4 v 4 |
| oxaliplatin | GSE227315 | **HCT116** (curated) | 0.5 µM, 24 h vs untreated | 3 v 3 |
| 5-fluorouracil | GSE343958 | **HCT116** (curated), 3D spheroids | 10 µM, 48 h vs DMSO | 4 v 4 |
| carboplatin | GSE185735 | SKOV3 (ovarian) | carboplatin vs control | 3 v 3 |
| erastin | GSE232034 | LNCaP | 5 µM vs DMSO | 3 v 3 |
| erastin | GSE232034 | PC3 | 5 µM vs DMSO | 3 v 3 |

Sample columns were read off GEO's own annotation, never inferred. GSE319411's
mapping, for example, comes from its per-sample `Library name: NC2_1` (control) and
`Library name: Bo2_1` (cisplatin) fields.

**The signatures reproduce known biology unprompted**, which is the best available
check that the ingest is correct:

- **oxaliplatin and 5-FU in HCT116** (p53 wild-type): top induced genes are CDKN1A
  (p21), MDM2, BTG2, FAS, TP53I3 — the canonical p53 damage response; repressed are
  RAD51, CCNE2, CENPI.
- **cisplatin in MDA-MB-231** (p53-mutant): no p21 induction, stress and
  inflammatory genes instead (EGR1, IL6, SPHK1, STC1). Consistent with the genotype.
- **erastin in LNCaP**: top induced is HMOX1, the canonical ferroptosis and
  oxidative-stress marker erastin is known to drive; repressed are the
  androgen-responsive KLK3 (PSA) and KLK2.

## The QC that reframes the result

Similarity between the drug signatures themselves, on 839 shared genes:

| comparison | cosine |
|---|---|
| two rapalogs, LINCS consensus (sirolimus vs everolimus) | **+0.87** |
| same cell line, two DNA-damaging drugs (oxaliplatin vs 5-FU in HCT116) | **+0.71** |
| mitomycin C vs doxorubicin, LINCS consensus | +0.42 |
| **same drug, two cell lines (erastin: LNCaP vs PC3)** | **+0.18** |
| two platinums, different cell lines (cisplatin MDA-MB-231 vs oxaliplatin HCT116) | −0.03 |

Erastin against erastin — same study, same dose, same platform, two cell lines —
reaches only 0.18, while two *different* drugs in one cell line reach 0.71. So the
response program is largely cell-line-specific. Two consequences:

1. **Cross-cell-line retrieval is weakly interpretable.** Almost every NPI × drug
   comparison available is cross-line, because only two pairs share a cell line.
2. **It explains why GSE95640 worked.** That retrieval used LINCS consensus
   signatures, which are medians over many cell lines and so average out exactly
   this cell-line component. Single-cell-line GEO signatures are real data but are
   not drop-in substitutes for consensus signatures in retrieval.

## Results, read with that caveat

14 human in-vitro cancer-line NPI signatures × 13 drug signatures, percentiles
against a 427-compound LINCS background. Full table in
`out/npi_drug_retrieval/npi_drug_similarity.tsv`.

**Consistent across three independent breast lines**, heat (GSE48398, 45 °C/30 min)
most resembles DNA-damaging chemotherapy:

| NPI | top matches (cosine, background percentile) |
|---|---|
| heat, MDA-MB-231 | oxaliplatin 0.24 (99th), doxorubicin 0.23 (97th), 5-FU 0.18 (95th) |
| heat, MCF7 | oxaliplatin 0.17 (98th), paclitaxel 0.17 (98th), doxorubicin 0.15 (95th) |
| heat, MDA-MB-468 | oxaliplatin 0.09 (92nd), temozolomide 0.08 (89th) |

**A mechanistically coherent match:** methionine deprivation (MCF7) most resembles
erastin (0.18, 87th percentile). Methionine feeds cysteine via transsulfuration and
so glutathione synthesis, which is the pool erastin depletes. If that holds, the two
duplicate a mechanism rather than adding two.

**The two same-cell-line comparisons, the only rigorous ones, show no mimicry:**
heat vs cisplatin in MDA-MB-231 = −0.04 (27th percentile); methionine deprivation
vs erastin in PC3 = −0.13 (15th). Both are the comparisons least affected by the
cell-line confound, and both are slightly negative — so the positive cross-line
similarities above should be treated as provisional.

**A caution against over-reading:** the curated tier-1 synergy is cystine
withdrawal + erastin, which act on the same transporter axis. Cystine deprivation's
similarity to erastin is only 0.03 (58th percentile). Transcriptomic similarity did
not detect a mechanism overlap that is certain on other grounds, which is one more
reason to treat these numbers as weak evidence.

## Deliberately excluded

- **Resistant-vs-parental comparisons** (GSE314768, GSE279527, GSE306865, GSE198042,
  GSE317900 and others), which scored well in the search but measure adaptation, not
  drug response.
- **triapine** GSE166419: clean design but one treated replicate per line, and
  Affymetrix transcript-cluster intensities with no probe map for GPL23126.
- **erastin in U-87** GSE236253: would add a curated cell line, but GEO does not map
  its featureCounts columns (`U87-3.bam`, …) to the WT and erastin samples. Not guessed.
- **lomustine**: no acute cell-line treatment series among 22 human expression GSEs.
- **cyclophosphamide**: none, for a mechanistic reason — a prodrug needing hepatic
  CYP activation, so cultured cells barely respond. Its curated interaction is an
  in vivo mouse result.

## What would make this stronger

1. **Build multi-cell-line consensus signatures for the platinums.** The search found
   many cisplatin series; averaging 3–5 independent lines would give a consensus
   comparable to LINCS and remove the cell-line confound that limits the results above.
2. **Pair NPI and drug in the same cell line.** Only two such pairs exist now. Heat
   or nutrient restriction profiled in HCT116, where oxaliplatin and 5-FU signatures
   already exist, would give matched comparisons.
3. **Do not read similarity as synergy** in either direction; this project measured
   that neither similarity nor orthogonality predicts it.

## Reproduce

```bash
python scripts/find_drug_datasets.py --missing-from-lincs --check-supplement
npi-pharma fetch-geo GSE319411 GSE227315 GSE343958 GSE185735 GSE232034 --suppl
python -c "from npi_pharma.ingest.drug_geo import build_all; from npi_pharma.store import save_signatures; \
  s,p=build_all('configs/drug_signatures_geo.yaml','data/raw'); print(p); save_signatures(s,'data/processed/signatures/drugs_geo.parquet')"
npi-pharma ingest-lincs --gctx ... --gene-space landmark \
  --perts metformin,sirolimus,everolimus,mitomycin-c,doxorubicin,paclitaxel,temozolomide \
  --out data/processed/signatures/lincs_curated.parquet
python scripts/npi_drug_retrieval.py --cancer-resource "<...>/Cancer Resource"
```
