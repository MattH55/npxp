# Mechanism-based NPI x drug hypotheses: real literature + real curated databases, not a model

Built 2026-10-07. Not a synergy-prediction model. Previous work in this
repo tested whether expression-signature similarity/composition predicts
real measured drug interaction, twice (`pipeline/docs/validation_drugcomb.md`
on 739,964 real DrugComb pairs, and a pathway-feature GBT in the separate
Cancer Resource repo on 2,924 real ALMANAC/DrugComb pairs), and found no
real signal either time. This does not attempt that again. It is a
curated lookup over real, independently-published mechanisms and real,
curated gene/drug databases, with every row traceable to specific
citations rather than a numeric score.

## The idea

A non-pharmacological intervention changes a cell's functional state.
Several of the modalities this project tracks have a real, independently
published mechanism describing *which* state, in *gene-level* terms:

| Modality | Real induced state | Citation |
|---|---|---|
| Mild hyperthermia (41-42.5C) | BRCA2 degradation, HR deficiency | Krawczyk et al., *PNAS* 2011, PMID 21555554 |
| Chronic hypoxia (<1% O2) | Decreased BRCA2/RAD51 synthesis, HR deficiency | Chan et al., *Cancer Research* 2008, PMID 18199558 |
| Tumor Treating Fields | BRCA1 downregulation, "BRCAness" | Giladi et al., *Cell Death & Disease* 2017, PMID 28358361 |
| Fasting / caloric restriction | IGF-1/PI3K/AKT/mTOR suppression | Nencioni et al., *Nature Reviews Cancer* 2018, PMID 30327499 |
| Exercise | Epinephrine/IL-6-dependent NK cell mobilisation | Pedersen et al., *Cell Metabolism* 2016, PMID 26895752 |

Three of these (hyperthermia, hypoxia, TTFields) independently converge on
the same real vulnerability -- a transient, non-mutational phenocopy of
BRCA-deficiency -- which is druggable via the same real, well-established
synthetic-lethal relationship BRCA-mutant cancers already exploit
(BRCA-PARP1, Bryant/Farmer, *Nature* 2005). That convergence is itself a
real, checkable finding, not something this pipeline derived -- three
separate wet-lab papers found it independently.

## What the pipeline does

For each modifier (`configs/npi_mechanisms.yaml`), by category:

- **induced_hr_deficiency** (hyperthermia, hypoxia, TTFields): look up the
  real induced lesion gene(s) (e.g. BRCA2) against
  [SynLethDB 3.0](https://zenodo.org/records/22843223) (37,943 real
  curated human synthetic-lethality pairs, CC-BY-4.0) for real SL
  partners, each with its own citation and evidence type. Then look up
  real approved anti-neoplastic drugs targeting each partner in
  [DGIdb 5.0](https://dgidb.org).
- **pathway_suppression** (fasting/CR): the literature's own recommended
  drug class (PI3K/AKT/mTOR inhibitors) is used directly; DGIdb is
  queried only to find real approved drugs in that class.
- **immune_mobilization** (exercise): real approved immune-checkpoint
  inhibitors (PDCD1/CD274/CTLA4/LAG3/HAVCR2/TIGIT), filtered to a direct
  inhibitor/antagonist/blocker claim.

### A real data-quality finding worth recording

DGIdb's own `immunotherapy` flag is broader than checkpoint blockade --
it also flags immunosuppressants used in a transplant context (sirolimus)
and various monoclonal antibodies (rituximab, cetuximab) that have
nothing to do with NK-cell/checkpoint biology. Likewise, querying a gene
like `CTLA4` or `PDCD1` directly returns many rows that are bare
biomarker/association claims (`interaction_type` is `NULL`) rather than
real drug-target relationships -- e.g. prednisone and cyclosporine show
up "associated with" CTLA4 in biomarker studies, not as CTLA4-targeting
drugs. `dgidb.py`'s `DIRECT_ACTION_TYPES` filter (require at least one
real source to report `inhibitor`/`antagonist`/`blocker`/etc., not only a
bare association) was added after checking this directly against the
real downloaded data -- without it, the checkpoint-gene lookups returned
mostly unrelated chemotherapies and immunosuppressants instead of the
real checkpoint inhibitors.

### Ranking

Rows are ranked by: (1) whether a *second*, independent real paper tested
this exact modifier with this exact drug target already (currently only
the Oncotarget heat+PARP-inhibitor+cisplatin follow-up, PMID 28427225,
and flagged only for PARP1 -- not broadcast to every other SL partner of
BRCA2 a drug happens to be found for, which an earlier version of this
script did incorrectly); then (2) SynLethDB's own evidence-type hierarchy
(CRISPR/CRISPRi > high/low-throughput screen or curated database >
computational prediction > text mining); then (3) how many independent
DGIdb source databases corroborate the drug-target claim.

## What this is not

Not a calibrated score. Not validated against any measured modifier x
drug outcome (see `pipeline/docs/validation_drugcomb.md` for why that
test, when run on the closest real proxy data this project has, found no
signal for the expression-based alternative). A row with many SL
partners and many corroborating drugs is not "more likely to work" in any
quantified sense -- it is a real, cited, mechanistically coherent
hypothesis, and the "direct experimental support" badge is the only
signal here that two independent studies actually agree on a specific
combination, which is a qualitatively different (and stronger) claim than
anything else on the page.

## Reproducing

```bash
cd pipeline
python scripts/fetch_mechanism_data.py        # real SynLethDB + DGIdb downloads
python scripts/build_mechanism_hypotheses.py --out out/mechanism_hypotheses.json
cp out/mechanism_hypotheses.json ../api/mechanism_hypotheses.json
```
