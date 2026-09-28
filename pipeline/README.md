# npi_pharma: NPI × drug interaction scoring for an individual transcriptome

This package scores how a non-pharmaceutical intervention (NPI) and a drug are
expected to interact **in one specific person**. It uses that person's
transcriptome plus precomputed intervention signatures. It is a
signature-composition engine in the style of CMap, comboSC and DIPx. It needs no
paired human NPI+drug labels.

> **Evidence tier:** every output is `tier_3_mechanism_only`, a plausibility
> ranking. It is not measured or calibrated synergy, not a PK prediction, and not
> clinical advice. Every report carries this disclaimer.
>
> That disclaimer now has a measured bound. Tested against DrugComb's 740k
> measured drug-pair outcomes, the composition scores do not predict interaction
> (see [docs/validation_drugcomb.md](docs/validation_drugcomb.md)). Read the
> scores as mechanism-only hypotheses, never as synergy estimates.

## Install

```bash
cd pipeline
pip install -e ".[test]"      # numpy, pandas, pyarrow, pyyaml (+ h5py, pytest)
python -m pytest              # 35 tests, ~2 s, no network
npi-pharma demo --out out/demo  # offline milestone path on SYNTHETIC fixtures
```

## Pipeline

```
configs/npi_catalog.yaml ─ ingest-npi ──► npis.parquet ─┐
LINCS Level 5 GCTX ─────── ingest-lincs ► drugs.parquet ─┼─ score / rank ─► report.json + ranked.tsv
patient expression ─────── encode-patient ► patient.npz ─┘
```

| module | role |
|---|---|
| `genes.py` | Normalises symbols and aligns every input onto one shared gene space. Raises `GeneOverlapError` below 200 genes instead of returning an empty space. |
| `ingest/lincs.py` | Reads only the requested columns of a Level 5 GCTX. Collapses each compound to a per-gene median, across all cell lines or a matched lineage. Resolves `rapamycin` ↔ `sirolimus`. `--random N --seed` draws a deterministic background set. |
| `ingest/geo.py`, `ingest/catalog.py` | Parses GEO series matrices, maps probes to genes, and splits pre/post from `!Sample_characteristics`. Validates catalog records: uncurated metadata is listed, never guessed. |
| `signatures/build.py` | `build_signature(pre, post, metadata)`: paired (or unpaired) moderated DE, then a robust z-score. Keeps the full vector plus the top 150 up/down genes. Consensus is refused across tissues or modalities. |
| `patient/encode.py` | Computes s_P as the patient's z-score against a healthy reference of the same tissue. Alternatively it takes the z against a same-batch patient cohort plus a within-study case-vs-healthy offset (`build-offset`, Cohen's d), flagged `health_via_external_offset`. Falls back to a cohort (flagged `cohort_relative`) or, only if forced, a within-sample z-score (flagged `no_reference`). Reads `.h5ad` with pseudobulk and an obs filter as the scTherapy-style hook. |
| `gene_sets.py` | Pathway z-activity `Σz/√k` (the same function for interventions and patients), plus ssGSEA. |
| `interact/score.py` | The v1 scores below. |
| `interact/rank.py` | Scores every NPI × drug pair and filters by tissue compatibility. Pairs that cannot be scored stay in the table with an `error:` flag. |

### Scores (`interact/score.py`)

s_P is the patient state to reverse; s_N and s_D are the NPI and drug
signatures. The combination is the sum of the two unit-normalised signatures,
s_N/|s_N| + s_D/|s_D|. Reversal is `reverse(P, s) = −cos(s_P, s)`, computed on
the patient's top 500 |s_P| "focus" genes.

| field | definition |
|---|---|
| `complementarity` | `gain = rev(combo) − max(rev(N), rev(D))`, divided by the remaining headroom `1 − max`. The raw value is kept as `complementarity_gain`. **Measured bound:** tested against 50,069 DrugComb drug-pair observations, this score has no detectable relation to measured synergy (residual Spearman ≤ 0.02 on four metrics) — see [docs/validation_drugcomb.md](docs/validation_drugcomb.md). |
| `orthogonality` | `1 − |cos(s_N, s_D)|` on focus genes. It counts toward the composite only when both agents reverse s_P. The ungated value is `orthogonality_raw`. |
| `monotherapy_correlation` | `cos(s_N, s_D)`. Used as a feature, following eLife 2020;9:e52707. |
| `pathway_joint` | DIPx-like. For each pathway, `need·(u_N+u_D)/2 − λ·|need|·conflict`, normalised by `Σ|need|`. |
| `mono_reversal` | Mean of `rev(N)` and `rev(D)`. |
| `composite` | `Σ wᵢ·scoreᵢ / Σ|wᵢ|`, times a confidence multiplier. |

All weights, shrink factors and thresholds are in
`src/npi_pharma/data/scoring.yaml`; override any of them with `--config`.

**Confidence multipliers:**

| condition | multiplier |
|---|---|
| tissue mismatch | ×0.5 |
| drug signature not from the patient's tissue (every LINCS cell line) | ×0.8 |
| NPI signature with n < 6 | ×0.6 |
| exploratory NPI | ×0.7 |
| uncurated NPI metadata | ×0.85 |
| cross-species comparison | ×0.8 |
| fewer than 1000 shared genes | ×0.8 |

**Flags (not multipliers):**
- `possible_redundancy`
- `pathway_antagonism`: NPI and drug push a pathway the patient needs in opposite directions.
- `xenobiotic_genes_moved_by_both`: raised only when the NPI induces xenobiotic-metabolism genes and the drug also moves them. It is a transcriptomic observation, not a PK claim.
- `SYNTHETIC_FIXTURE`

Attributions come in two forms. Per-gene terms of `rev(combo)` are listed as
`genes_driving_reversal` and `genes_opposing_reversal`. Per-pathway joint terms
are listed under `pathways`.

## First milestone on real data

The milestone has been run on GSE95640 (LCD, adipose RNA-seq, n = 191 pairs) against
the LINCS Phase II L1000 signatures. See [docs/milestone_GSE95640.md](docs/milestone_GSE95640.md)
for the inputs, the curation checks, the result, and what the result does not show.

The follow-up in the same doc adds four things:
- a health-referenced s_P (GSE244118 obese-vs-lean offset)
- held-out signatures (`ingest-npi --exclude-samples`)
- a 191-subject sweep (`scripts/cohort_sweep.py`)
- a regression check against the paper's L1000CDS2 hits (`scripts/retrieval_check.py`): AUC 0.75, p = 0.0005
To reproduce:

```bash
# 1. Drugs: LINCS Phase II Level 5 from GEO GSE70138 (5.4 GB download; unpacked GCTX 5.8 GB)
npi-pharma fetch-lincs --gse GSE70138 --raw-dir data/raw      # prints the exact ingest-lincs command
npi-pharma ingest-lincs --gctx data/raw/GSE70138/<...>.gctx --sig-info ... --gene-info ... \
  --drug-set metabolic --random 20 --seed 0 --out data/processed/signatures/drugs.parquet

# 2. NPI. GSE95640 is RNA-seq: its counts are in the supplementary files, and
#    --gene-info writes the Ensembl -> symbol map from NCBI gene_info.
npi-pharma fetch-geo GSE95640 --suppl --gene-info --raw-dir data/raw
python scripts/curate_gse95640.py        # subject pairing + timepoint checks -> samples_curated.tsv
npi-pharma ingest-npi --ids LCD_adipose_GSE95640
#    Microarray series (e.g. GSE77962) instead use the series matrix + GPL probe map:
#    npi-pharma inspect-geo data/raw/GSE77962/GSE77962_series_matrix.txt.gz, then fill inputs.pairing

# 3. Patient: one adipose baseline, against healthy adipose references of the same platform
npi-pharma encode-patient --input patient.tsv --sample GSMxxxx --tissue adipose --reference healthy.tsv --out patient.npz

# 4. Score + ranked table + JSON explanation of which pathways drove the top pair
npi-pharma score --patient patient.npz --npi LCD_adipose_GSE95640 --drug metformin --out report.json --tsv ranked.tsv
npi-pharma rank  --patient patient.npz --npi-class diet --drug-set metabolic --top 25 --out rank.tsv --json rank.json
```

`fetch-*` need HTTPS access to `ftp.ncbi.nlm.nih.gov`.

Catalog `inputs` accept three forms:
- `series_matrix` + `probe_map` (microarray)
- `counts` + `id_map` + `samples` (RNA-seq raw counts, converted to log2(CPM+1) with a low-expression filter)
- `expression` + `samples` (a genes × samples matrix that is already normalised)

For real Hallmark/KEGG analyses, pass `--gene-sets h.all.v2024.1.Hs.symbols.gmt`.
The bundled `curated_core.gmt` holds 12 compact, hand-curated core sets
(mTORC1, insulin, inflammation, OXPHOS, xenobiotic, FAO, lipogenesis,
cholesterol, hypoxia, UPR/ISR, AMPK, E2F). They are not MSigDB.

## Fixtures

`npi_pharma/fixtures.py` generates a deterministic **synthetic** world: 50
drug signatures and 5 NPI signatures built from latent pathway programs over
the bundled gene sets. The CR and exercise arms are built through the real
`build_signature` from mouse-case pre/post matrices. The world also includes an
adipose patient and a healthy reference. Every id starts with `FIXTURE_`, every
record has `provenance: synthetic_fixture`, and every score carries the
`SYNTHETIC_FIXTURE` flag. Unit tests never download LINCS.

## Not yet built (by design)

The brief stops at the milestone. This package does not yet include:
- scRNA-seq clone-wise DE, beyond the pseudobulk hook
- NPI dose–response modelling
- PK simulation
- the v2 trainable head (drug–drug pretraining followed by an NPI encoder)

v2 must not report cell-line drug-synergy percentages as if they were NPI × drug
measurements.
