# References

Papers and resources the `npi_pharma` pipeline is built around. Accessions used
by signature records live in `configs/npi_catalog.yaml`.

## Drug × individual transcriptome (frameworks this extends)
- PDSP: cell-line pretraining, then fine-tuning on patient expression. https://academic.oup.com/bioinformatics/article/40/5/btae134/7667298 · https://github.com/hikuru/PDSP
- scTherapy: malignant/normal split, then per-clone drug prediction. https://www.nature.com/articles/s41467-024-52980-5 · https://github.com/kris-nader/scTherapy
- comboSC: patient scRNA-seq × CMap combinations. https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-023-01256-6 · https://github.com/bm2-lab/comboSC
- DIPx: pathway-level synergy explanations (the model for `pathway_joint`). https://elifesciences.org/articles/100071
- AML scRNA-seq + XGBoost combinations. https://aacrjournals.org/cancerres/article/85/14/2753/763446

## Signature math
- CMap / L1000 (Subramanian et al., Cell 2017). https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5990023/
- A combination's transcriptome is not the sum of its monotherapies; correlated monotherapy programs go with synergy. https://elifesciences.org/articles/52707

## LINCS L1000 data access
- GEO Phase I GSE92742: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742
- GEO Phase II GSE70138: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138
- clue.io: https://clue.io · data dashboard https://clue.io/releases/data-dashboard
- Where to get L1000: https://lincsproject.org/LINCS/tools/workflows/find-the-best-place-to-obtain-the-lincs-l1000-data
- Extraction helper: https://github.com/BenderGroup/LINCS-Extraction
- Optional signature search (not a hard dependency): https://maayanlab.cloud/sigcom-lincs/ · https://maayanlab.cloud/L1000CDS2/

## NPI × pharma, or NPI signatures
- CR + rapamycin joint RNA-seq (regression reference for combo ≠ max(mono)). https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12078523/
- LCD adipose signature (GSE95640) → L1000CDS2 drug mimics (regression reference for NPI → drug retrieval). https://link.springer.com/article/10.1186/s12967-025-07424-z
- LCD vs VLCD adipose time course. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77962
- Fucoxanthin + PI3K inhibitor. https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0239551
- Human muscle diet + exercise, individual response and G×L eQTLs. https://www.cell.com/cell-genomics/fulltext/S2666-979X(25)00207-1
- Exercise vs obesity, tissue-specific single-cell programs. https://www.sciencedirect.com/science/article/pii/S1550413122003941
- CR vs exercise vs both (mouse RNA-seq). https://www.mdpi.com/2072-6643/15/4/1047 · mammary CR vs exercise https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2798590
- Baseline expression predicts response to a meditation-class NPI. https://www.sciencedirect.com/science/article/abs/pii/S0889159124005269
- PGx + nutrigenomics framing. https://www.mdpi.com/2075-4426/14/12/1121

## Synergy prediction from gene expression (what the field has, and how it evaluates)

Models that predict drug-pair synergy from **baseline (untreated) cell-line
expression** plus drug chemical features. None of them takes a post-treatment
signature as the input, and none takes an NPI as an agent.

- DeepSynergy (Bioinformatics 2018): the original DNN on chemical descriptors +
  cell-line expression. https://academic.oup.com/bioinformatics/article/34/9/1538/4747884
- MatchMaker: two per-drug subnetworks on DrugComb; reports Pearson 0.79,
  Spearman 0.74. https://www.biorxiv.org/content/10.1101/2020.05.24.113241v1.full
- DRSPRING (GCN-based). https://www.sciencedirect.com/science/article/pii/S0010482524005201
- PerturbSynX: the nearest relative to what this project wants, in that it uses
  perturbation rather than only baseline data.
  https://www.sciencedirect.com/science/article/abs/pii/S1476927125003974
- 2024 mini-review: best models do well on *known* drugs and cell lines (AUROC to
  0.98) while "scenarios involving new drugs or cell lines still fall short";
  leave-drug-out is much harder than leave-pair-out. https://arxiv.org/html/2404.02484
- Systematic evaluation (PLOS Comput Biol 2022).
  https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010200

**The decisive one for this project.** SynVerse (Brief Bioinform 2025) evaluated 16
synergy models over 8 drug/cell-line feature types, 5 preprocessing schemes and 2
encoders, under four splitting strategies, with module, feature-shuffling and
network-based ablations. **None outperformed a naive one-hot-encoding baseline**, and
models with *shuffled* drug and cell-line features performed comparably to those with
the real ones -- so performance was not driven by biologically informative features.
https://academic.oup.com/bib/article/26/6/bbaf676/8407512 ·
https://pmc.ncbi.nlm.nih.gov/articles/PMC12753315

That is this project's own measurement arrived at independently: signature
composition carries no detectable synergy information once drug and cell-line main
effects are removed (docs/validation_drugcomb.md, <= 0.02 residual Spearman on
739,964 observations against a -0.30 positive control). The negative result here is
the field's result under proper ablation, not a failure of this implementation.

## Four families, and which could take an NPI as an agent

I earlier said there is "no transfer" from synergy modelling to NPI x drug. That was
too broad: it was true of one family and stated as if it covered the field. There are
four, and they fail or transfer for different reasons.

**1. Baseline expression + chemical structure** (DeepSynergy, MatchMaker, DRSPRING;
above). Cannot take an NPI at all -- the drug encoder consumes SMILES and hyperthermia
has no structure -- and SynVerse's shuffling ablation shows the accuracy rests on agent
identity seen in training. Measured here independently: cell-line expression is
interchangeable with randomly reassigned expression (`validation_drugcomb.md`).

**2. Perturbation-composition models.** CPA decomposes expression into additive basal,
perturbation and covariate embeddings and predicts **unseen perturbation
combinations** by swapping embeddings, including combinatorial genetic interactions;
GEARS predicts outcomes of novel multi-gene perturbations from a gene-gene graph.
*This family is structurally right for NPIs*, and that is the part I was wrong to
dismiss: a perturbagen is identified by its measured effect, not by a chemical
structure, and GEARS's perturbagens are genes, which have no SMILES either. An NPI is
exactly that kind of agent.
- CPA: https://www.biorxiv.org/content/10.1101/2021.04.14.439903v1 · https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10258562/
- GEARS: https://www.nature.com/articles/s41587-023-01905-6 · https://github.com/snap-stanford/GEARS
- MultiCPA: https://www.biorxiv.org/content/10.1101/2022.07.08.499049.full.pdf

  **But this family has its own SynVerse, and it is just as blunt.** "Deep-learning-based
  gene perturbation effect prediction does not yet outperform simple linear baselines"
  (Nature Methods 2025) -- no deep model consistently beat the mean prediction or a
  linear model. PerturBench found simple baselines matching or beating GEARS and
  Geneformer, and for **unseen two-gene perturbations the matching-mean baseline beat
  GEARS by 11% on Pearson delta**. Unseen combinations are precisely the regime an NPI x
  drug application would live in.
  - https://www.nature.com/articles/s41592-025-02772-6
  - PerturBench: https://arxiv.org/pdf/2408.10609
  - Systema: https://www.nature.com/articles/s41587-025-02777-8

**3. Monotherapy-signature correlation.** The eLife result already cited above is a
*predictor*, not only an observation: drug synergy tends to follow when two drugs'
individual transcriptional responses are correlated, and a simple algorithm on that
basis reached **AUROC 0.77** on an independent dataset.
https://elifesciences.org/articles/52707

  This transfers directly and needs nothing an NPI lacks -- two monotherapy signatures,
  both of which this project has. **It has also already been tested here and it did not
  hold**: signature similarity and complementarity carry <= 0.02 residual Spearman
  against 739,964 DrugComb pairs once drug and cell-line main effects are removed,
  against a -0.30 positive control (`validation_drugcomb.md`). The likeliest
  reconciliation is that the eLife evaluation is a smaller independent set without
  main-effect ablation, which is the same gap SynVerse found across the other family.
  Worth re-reading before any further work on this axis, because if it is right the
  disagreement is informative and if it is not, the 0.77 is a main-effect artefact.

**4. Pathway-activation and mechanistic models.** DIPx and the pathway-activation
model (both cited above) predict synergy from pathway activity rather than from drug
identity. NPIs map onto pathways cleanly (heat -> HSF1/HSP90, fasting -> AMPK/mTOR,
cystine deprivation -> SLC7A11/GPX4), so this family can take an NPI without
modification. It predicts direction and mechanism rather than a calibrated synergy
score, which is the honest ceiling with no labels.
Also relevant: disease-specific drug response profiles plus single-cell signatures,
https://elifesciences.org/reviewed-preprints/102442

### Summary

| family | can take an NPI? | why it is still blocked |
|---|---|---|
| baseline expression + structure | no | no SMILES; identity-driven accuracy |
| perturbation composition (CPA/GEARS) | **yes, structurally** | the family does not beat mean/linear baselines on unseen combinations |
| monotherapy-signature correlation | **yes, directly** | measured here at <= 0.02 residual after ablation |
| pathway / mechanistic | yes | gives direction, not a calibrated score |

So the accurate statement is not "no transfer". It is that two of the four families
*can* accept an NPI as an agent, and both are blocked by evidence rather than by
plumbing -- one by its own benchmarking literature, one by this project's own
measurement.

## Later (v2 trainable head)
- PAIRWISE (transcriptome + drug graphs + targets). https://github.com/Mew233/pairwise
- JointSyn (dual-view, sample-specific synergy; Bioinformatics 2024).
- MARSY is already named as the Stage 3 base model in `../../PREDICTION_METHODOLOGY.md`.
