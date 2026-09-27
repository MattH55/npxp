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

## Later (v2 trainable head)
- PAIRWISE (transcriptome + drug graphs + targets). https://github.com/Mew233/pairwise
- JointSyn (dual-view, sample-specific synergy; Bioinformatics 2024).
- MARSY is already named as the Stage 3 base model in `../../PREDICTION_METHODOLOGY.md`.
