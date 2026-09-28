# Search log: NPI single-agent survival per cell line

2026-09-28. Following the precedent of `Cancer Resource/data/modifier_pairs/PAIR_SEARCH_LOG.md`,
this records a search that mostly failed, so the next attempt does not repeat it.

**Target.** `configs/npi_monotherapy.yaml` needs, for each NPI and each curated
cell line, a published surviving fraction or viability loss for the **NPI alone**.
That is the input the efficacy predictor needs
([efficacy_from_monotherapy.md](efficacy_from_monotherapy.md)) and the quantity
the literature identifies as the strongest predictor of a combination's outcome.
The thermal enhancement ratios already curated are a different quantity
(enhancement of a drug) and cannot substitute.

**Outcome: 3 usable entries over 2 cell lines, all from one paper** (Helderman
et al. 2020, found before this search). Five systematic strategies added nothing.

## What was searched

`scripts/find_npi_monotherapy.py` queries PMC open-access full text through NCBI
E-utilities, splits articles into sentences, and keeps those naming a curated cell
line together with a number and a survival/viability word. It prints candidates
for a human to judge and never writes the YAML itself.

| pass | query focus | papers scanned | candidate sentences | usable |
|---|---|---|---|---|
| 1 | hyperthermia + surviving fraction / clonogenic | 146 | 51 | 0 |
| 2 | water-bath hyperthermia, HIPEC, 41–43 °C + clonogenic | 116 | 12 | 0 |
| 3 | glucose restriction / deprivation + viability | 40 | 3 | 0 |
| 4 | cystine, cysteine, methionine deprivation + viability | 40 | 0 | 0 |
| 5 | PubMed **abstracts** (all years, not only open access) | 28 | 2 | 0 |

## Why it failed, specifically

1. **The modality collides.** "Hyperthermia" in recent open-access literature is
   overwhelmingly nanoparticle photothermal therapy — gold nanorods, magnetic
   hyperthermia, photodynamic therapy. Those papers report viability numbers, but
   for localised nanoparticle heating with a carrier, not the water-bath
   whole-culture heating that HIPEC and thermoradiotherapy studies use. The
   numbers are not transferable, and they dominate every keyword search.
2. **The classic literature is not in the open-access full text.** The thermal
   biology that established these survival curves (the work behind the CEM43
   thermal-dose concept) is indexed in PubMed but its full text is not in the PMC
   open-access subset, so sentence-level search cannot reach it. Pass 5 searched
   abstracts to get around this; abstracts state conclusions, not per-line
   survival percentages.
3. **The numbers live in figures.** Survival is reported as a curve or a heat map.
   The Helderman extraction worked only because two sentences happened to state
   values in words ("only 19% of the cells surviving", "a 30% reduction"). The same
   paper's Figures 2–5 are colour gradients with no printed values, which its own
   curation note in `data/thermal_ter/ter_corpus.json` already recorded.
4. **A local copy was not what it claimed.** `Cancer Resource/data/papers/hipec_2020/supplement.pdf`
   is an HTML error page saved with a .pdf extension, not the supplement.

## A proxy was tested and is not good enough

DepMap CRISPR dependency is a tempting substitute: a line that dies without
SLC7A11 is a line that dies without cystine import, which is what cystine
withdrawal does. `scripts/test_dependency_proxy.py` calibrates the idea against
cases where the link is textbook — a gene whose knockout phenocopies its own
inhibitor — using the best of PRISM, GDSC1 and GDSC2 (expected sign positive: less
dependent, higher AUC, more resistant):

| calibration pair | best Spearman |
|---|---|
| EGFR × gefitinib | **+0.47** |
| MDM2 × nutlin-3 | +0.43 |
| EGFR × erlotinib | +0.46 |
| BCL2L1 × navitoclax | +0.30 |
| BRAF × dabrafenib | +0.22 |
| CDK4 × palbociclib | +0.22 |
| ABL1 × imatinib | +0.07 |
| FLT3 × quizartinib | +0.01 |

| NPI proxy pair | best Spearman |
|---|---|
| SLC7A11 × erastin | +0.09 |
| MTOR × everolimus | +0.09 |
| MTOR × sirolimus | +0.08 |
| SLC7A11 × sulfasalazine | −0.01 |

Median of the best source: **+0.217 for calibration, +0.084 for the NPI proxies.**
So the approach can work — it reaches +0.2 to +0.47 when a dependency genuinely
drives a drug's effect — but the transporter proxies sit at the level of its
failures. Gene knockout is not nutrient withdrawal: it is permanent rather than
transient, complete rather than partial, and leaves time for compensation.
Not used.

**Useful by-product:** this is a calibration for the whole project. A per-cell-line
functional measurement predicting single-agent sensitivity tops out near +0.2 to
+0.47 even in the most favourable case available. Any feature claiming more than
that for a harder target deserves scrutiny.

## What would actually unblock this

1. **Full-text access to specific papers, then figure and table reading.** This is
   a library-access problem, not a search problem. The searches above identify the
   journals that publish it — *International Journal of Hyperthermia*, *Cells*,
   *Cancers* — and thermoradiotherapy work in cervical lines (HeLa, CaSki, SiHa)
   appeared in pass 5 (PMID 39837264) with clonogenic survival at 42 °C/1 h, but
   behind full text this session cannot read.
2. **Measure it.** For a lab holding these lines, heat-alone clonogenic survival at
   41/42/43 °C for 60 min is a standard experiment, and it is the single
   highest-value input for this project: 11 cell lines × 3 temperatures fills the
   whole thermal arm.
3. **Pick NPIs with systematic screens.** If a published screen reports a
   metabolic-stress sensitivity across many CCLE lines, that populates the NPI
   side at scale in one ingest, the way PRISM does for drugs.

## Reproduce

```bash
python scripts/find_npi_monotherapy.py --npi hyperthermia --max-papers 50
python scripts/find_npi_monotherapy.py --npi glucose_restriction --require-alone
python scripts/find_npi_monotherapy.py --npi amino_acid_restriction --require-alone
python scripts/test_dependency_proxy.py      # needs CRISPRGeneEffect.csv (figshare 27993248)
```
