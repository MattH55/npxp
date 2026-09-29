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

## Next

A re-run with the RNA-seq path enabled is the test of whether the platinums can reach
the series counts the validation says they need. The same two rules decide the
outcome: series count, then cross-series agreement. Any consensus below about 0.1
cross-series agreement should be treated as unusable regardless of how many series
went into it.

## Reproduce

```bash
python scripts/build_drug_consensus.py                      # the drugs LINCS lacks
python scripts/build_drug_consensus.py --drugs doxorubicin paclitaxel temozolomide \
  metformin "mitomycin C" etoposide gemcitabine vorinostat tamoxifen bortezomib \
  --out data/processed/signatures/drugs_consensus_validation.parquet \
  --report out/drug_consensus_validation                    # validation against LINCS
npi-pharma ingest-lincs --gctx ... --gene-space landmark --perts @perts.txt \
  --out data/processed/signatures/lincs_all.parquet         # all 1,763 LINCS compounds
```
