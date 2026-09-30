# Multi-cell-line consensus drug signatures, and how far they can be trusted

Run 2026-09-28/29. Code: `scripts/build_drug_consensus.py`,
`probe_map_from_platform_table` in `src/npi_pharma/ingest/fetch.py`.

A single-cell-line drug signature is dominated by the cell line, not the drug: in
this project's own measurement, two *different* DNA-damaging drugs in one cell line
agree at cosine 0.71 while the *same* drug in two cell lines agrees at 0.18
([npi_drug_retrieval.md](npi_drug_retrieval.md)). LINCS avoids that by taking a
median over many cell lines, which is why the one validated retrieval result used
LINCS consensus signatures. So every drug needs a consensus.

## Coverage

| source | drugs | how the consensus is formed |
|---|---|---|
| LINCS Phase II (GSE70138) | **1,763** | already a median over a median of 7 cell lines per compound |
| GEO series, assembled here | 10 (of 20 attempted) | per-gene median over independent series |
| LINCS Phase I (GSE92742) | ~19,800 | **not usable here**: 21.3 GB against 14 GB of free disk |

LINCS Phase II lacks every platinum agent, 5-fluorouracil, erastin, triapine,
lomustine, cyclophosphamide, methotrexate, pemetrexed, topotecan, the vinca
alkaloids, melphalan, carmustine, bleomycin, actinomycin D and 5-azacytidine — 22 of
54 widely used drugs checked. Those are the ones a GEO consensus has to cover.

## Validation: does the method recover a trusted consensus?

Ten drugs that LINCS *does* carry were run through the same GEO pipeline, so the
result can be compared against a consensus built from far more data.

| drug | GEO series used | cross-series agreement | **cosine vs LINCS consensus** |
|---|---|---|---|
| vorinostat | 8 | 0.119 | **0.652** |
| bortezomib | 4 | 0.334 | **0.554** |
| doxorubicin | 5 | 0.111 | **0.522** |
| etoposide | 3 | 0.098 | **0.396** |
| gemcitabine | 3 | 0.070 | **0.381** |
| tamoxifen | 3 | −0.091 | 0.154 |
| temozolomide | 2 | −0.012 | 0.015 |
| metformin | 3 | 0.014 | −0.065 |

Median 0.388; 5 of 8 reach 0.3 or better. For scale, two independent LINCS
consensus signatures of the same mechanism (sirolimus vs everolimus) agree at 0.87,
and two single-cell-line signatures of the same drug agree at 0.18. So a GEO
consensus built from enough series lands well above single-series noise and well
below a same-mechanism ceiling.

### Two usable rules come out of this

1. **More series, better consensus**: series count correlates 0.74 with agreement
   against LINCS. Two series is not enough; five to eight works.
2. **Cross-series agreement predicts quality without needing ground truth**: it
   correlates 0.70 with the LINCS comparison. That makes it a trust metric for drugs
   where no LINCS signature exists to check against.

Rule 2 is confirmed by its failures as well as its successes: the three drugs that
failed validation (metformin, temozolomide, tamoxifen) are exactly the three whose
cross-series agreement was at or below zero.

## Applying the trust rule to the drugs LINCS lacks

| consensus | series | cross-series agreement | verdict |
|---|---|---|---|
| erastin | 2 | 0.683 | inflated — both arms come from one study (LNCaP and PC3 of GSE232034), not independent |
| actinomycin D | 3 | 0.347 | usable |
| 5-fluorouracil | 2 | 0.219 | borderline |
| cyclophosphamide | 2 | 0.120 | weak |
| melphalan | 2 | 0.107 | weak |
| cisplatin | 2 | 0.016 | **not usable** |
| topotecan | 2 | 0.016 | **not usable** |
| 5-azacytidine | 3 | 0.008 | **not usable** |
| carboplatin | 2 | −0.040 | **not usable** |
| methotrexate | 2 | −0.049 | **not usable** |

So the platinums — the drugs carrying most of this project's tier-1 interaction
evidence — do **not** yet have a trustworthy consensus. They are built from only two
series each, which the validation above shows is too few.

## Why so few series, and what was fixed

Of 194 candidate series examined in the first pass, 164 were unusable. Four
separate causes, all now addressed:

1. **No probe map.** GEO's `.annot.gz` exists for major platforms (GPL570, GPL96,
   GPL6244, GPL10558, GPL6480) but not for newer Affymetrix and Agilent ones
   (GPL15207, GPL21185). Without it, rows stay as probe ids and a consensus shares no
   genes at all — silently. `probe_map_from_platform_table` falls back to the
   submitted platform table, parsing either a symbol column or Affymetrix's
   `gene_assignment` format ("NM_x // SYMBOL // …").
2. **RNA-seq series.** Their series matrix carries no expression table; counts sit in
   supplementary files under arbitrary column names. `match_columns` maps samples to
   columns using labels GEO itself publishes — the title, the `description` field
   (often an explicit "Library name: NC2_1"), and the supplementary-file basename —
   and refuses anything ambiguous rather than guessing.
3. **The CSV field-size limit.** GEO tables embed whole GO/KEGG blobs in single
   fields, past the 128 KiB default, which aborted the read outright.
4. **Wrong identifier column.** Tables commonly carry both an Ensembl id and a symbol
   column; the column is now chosen by content (which best matches known gene
   symbols), since only symbols compare across series.

A fifth problem was a correctness issue rather than a coverage one: arm detection was
selecting **resistant-versus-parental** designs, which measure adaptation rather than
drug response. Those are now rejected by field name and value, as are combination arms.

## The re-run has been done. The answer is no.

The v2 run enabled the RNA-seq ingest path that v1 discarded 164 of 194 candidate
series for. It produced **12 consensus signatures from 41 series** across the 20 drugs
LINCS Phase II lacks; the other 8 drugs found fewer than 2 usable series at all.

| drug | series | cross-series agreement | band |
|---|---|---|---|
| 5-azacytidine | 2 | 0.191 | low |
| 5-fluorouracil | **7** | 0.176 | low |
| erastin | **7** | 0.140 | low |
| vinblastine | 2 | 0.094 | very low |
| oxaliplatin | 3 | 0.087 | very low |
| carboplatin | 4 | 0.052 | very low |
| bleomycin | 3 | 0.033 | very low |
| topotecan | 2 | 0.016 | very low |
| cisplatin | 3 | 0.009 | very low |
| methotrexate | 2 | −0.017 | very low |
| melphalan | 3 | −0.027 | very low |
| vincristine | 2 | **−0.161** | very low |

**Not one drug reaches 0.30**, the threshold this document's own validation set for
"medium". Three reach 0.10. Vincristine is the sharpest statement of the problem: two
independent series of the same drug *anti-correlate*.

The re-run did what it was built for on **coverage** — 5-FU and erastin went from too
few series to 7 each, clearing the 5–8 series bar. It did nothing for **agreement**,
which is the quantity that actually predicts trustworthiness (r = 0.70 against the
LINCS validation). More series of a weak signal is still a weak signal.

For the platinums specifically, the drugs this whole exercise targeted: cisplatin
0.009, carboplatin 0.052, oxaliplatin 0.087. All below even the "low" floor. **A GEO
consensus is not a route to usable platinum signatures**, and anything built on one
rests on essentially nothing. Per the instruction to flag rather than halt, they are
computed, stored and banded — but the band is the result.

### What would actually work

Not more GEO series. Phase I was the obvious candidate and has now been scoped from
its metadata (`scripts/scope_lincs_phase1.py`), which changes the picture in both
directions.

**Correction: Phase I does not carry the platinums either.** An earlier note in this
project said it did. It carries **cisplatin only, in 4 cell lines** — fewer than the 7
that make a Phase II consensus trustworthy. **Carboplatin and oxaliplatin are absent
from Phase I entirely**, as are bleomycin, carmustine, lomustine and triapine. So the
platinum arm cannot be rescued from LINCS at all, and the conclusion above stands and
hardens: it needs a dedicated experiment or it does not get done.

**But Phase I is worth having for other reasons.** Ten drugs the project needs today
would reach "high" reliability, with cell-line counts far beyond anything GEO gave:

| drug | Phase I cell lines | today |
|---|---|---|
| erastin | **51** | GEO consensus, 7 series, agreement 0.140 ("low") |
| idarubicin | **51** | nothing |
| pemetrexed | **50** | nothing |
| vinblastine | 15 | GEO consensus, 2 series, 0.094 |
| vincristine | 14 | GEO consensus, 2 series, **−0.161** |
| vinorelbine | 14 | nothing |
| topotecan | 13 | GEO consensus, 2 series, 0.016 |
| dactinomycin | 12 | GEO consensus (actinomycin D) |
| methotrexate | 8 | GEO consensus, 2 series, −0.017 |
| cyclophosphamide | 8 | nothing |

The erastin row is the one that matters most. Cystine deprivation retrieving erastin
first of 1,803 signatures is this project's single best result, and it currently rests
on a consensus banded "low". Phase I would replace that with a median over 51 cell
lines — turning the best result into one whose input is trustworthy on the same terms
as the rest of the LINCS panel. Vincristine is the other notable row: its GEO
consensus is *negatively* self-consistent, and Phase I has 14 lines.

**It does not fit in this container.** The Level 5 matrix is 19.9 GB gzipped against
13.4 GB free, and there is no partial fetch: a `.gctx` is HDF5, which needs random
access, and a gzip stream cannot be seeked, so no subset of the remote file can be
read. Reading it also needs the decompressed HDF5, several times the archive, so
19.9 GB is a floor on the space required rather than the total.

Two routes that would work, neither of them a download from here:

* Run the ingest on a machine with ~80 GB free and commit the resulting parquet, which
  is small — the existing `npi-pharma ingest-lincs` already does the subsetting, and
  only ~500 signatures over 978 landmark genes are actually wanted.
* Use the CLUE API (`api.clue.io`), which serves Phase I signatures individually and
  would need a few MB rather than 20 GB. It requires a free registered key, which is
  the user's to obtain; nothing here can be done without it.

## Three defects found while reading these numbers

Each one was a silent failure that produced a plausible-looking result, which is why
they are recorded here rather than only in the commit log.

**A parse failure that arrived as zeros.** GSE304295 ships its TPM table under a
European locale (semicolon-separated, `2,74724752907299`). `to_numeric` made every
value NaN, and `groupby().sum()` of an all-NaN group is `0.0`, not NaN — so the
failure became a table of zeros. It passed the row count, the symbol check and every
other gate, gave a completely flat signature, and standardising that later produced
0/0, which turned oxaliplatin's median over three series into **NaN**. Downstream,
`reliability_for` bands a NaN exactly as it bands a measured 0.05, so a broken input
was about to be reported as a weak result. Fixed at all three points: comma decimals
are retried, `is_degenerate` rejects a constant table whatever the cause, and the QC
median ignores non-finite pairs and reports `n_pairs_computable`.

**An identifier that was not a column.** With the TPM file rejected, the reader fell
through to the same series' raw-counts file — which has no header over its gene
column, so pandas reads the gene names as the *index*. The identifier is chosen by
content rather than by name, which was right, but only columns were scored; it picked
a sample column and the symbol gate correctly rejected the series at 0% matched. The
index is now scored on the same footing. GSE304295 then builds properly (27,296
genes, 6 vs 6), which is why oxaliplatin has three series above rather than two.

**A cross-panel comparison that inflated everything.** `rank_npi_drug_all.py` reported
each row's percentile against the LINCS panel, including for GEO rows. The panels are
not on one scale:

| panel | mean \|similarity\| |
|---|---|
| LINCS consensus | 0.060 |
| GEO consensus | 0.109 |
| GEO consensus v2 | 0.099 |
| GEO factorial main effect | 0.115 |
| GEO single series | 0.118 |

A LINCS consensus is a median over ~7 cell lines, which shrinks the cell-line-specific
component; a GEO single- or few-series signature keeps it, and the NPI signatures are
themselves single-cell-line, so they share it. For the three amino-acid-deprivation
NPIs the gap is tenfold (0.008 against 0.07–0.12) and has nothing to do with the
drugs. Every GEO row's percentile against LINCS was therefore inflated by
construction. `rank_in_panel`, `percentile_in_panel` and `z_in_panel` are now the
figures to read, and the cross-panel column carries a warning in `summary.json`.

## One genuine positive control came out of it

**Cystine deprivation in MCF7 retrieves erastin first, out of all 1,803 drug
signatures** (similarity 0.338; rank 1 of 12 within its own panel, z = 2.83).

That is mechanistically exact and was not planted. Erastin induces ferroptosis by
inhibiting SLC7A11, the cystine/glutamate antiporter; cystine deprivation removes the
substrate that same transporter imports. Two routes to one mechanism, and the
retrieval put them together without being told. Methionine deprivation ranks erastin
third, which is also coherent — methionine feeds cysteine through transsulfuration,
so restricting it depletes the same axis more weakly.

Note the tension worth keeping in view: **this best-in-project result rests on a
signature banded "low"** (erastin's cross-series agreement is 0.140). The band is a
conservative proxy for trustworthiness, not a verdict on usefulness. It should not be
read as licence to trust the other low-agreement consensuses — one apt result does
not validate a panel — but it does show the band can be pessimistic about a signature
that is capturing real biology.

Before the within-panel correction above, this result looked like "32.5 SD above the
LINCS panel mean". That number was an artefact of the panel-scale confound. The
finding survived the correction; the overstatement did not.

## Reproduce

```bash
python scripts/build_drug_consensus.py                      # the drugs LINCS lacks
python scripts/build_drug_consensus.py --drugs doxorubicin paclitaxel temozolomide \
  metformin "mitomycin C" etoposide gemcitabine vorinostat tamoxifen bortezomib \
  --out data/processed/signatures/drugs_consensus_validation.parquet \
  --report out/drug_consensus_validation                    # validation against LINCS
python scripts/build_drug_consensus_all.py                  # the v2 re-run, one process per drug
python scripts/scope_lincs_phase1.py                        # what Phase I would add, from 12 MB of metadata
npi-pharma ingest-lincs --gctx ... --gene-space landmark --perts @perts.txt \
  --out data/processed/signatures/lincs_all.parquet         # all 1,763 LINCS compounds
```
