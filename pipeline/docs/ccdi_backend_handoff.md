# Hand-off: proposed backend table changes from the CCDI integration

These files are produced by `npi-pharma build-ccdi` from OpenPedCan v15
(CC-BY-4.0). The changes they imply belong in the backend workspace
(`synlethality/`, exported by `export_static.py`). None of them was made by
hand-editing `api/`. See `ccdi_integration.md` for methods and caveats.

## 1. New pediatric `cell_line` rows (candidates)

Source: `docs/ccdi/model_concordance_openpedcan_v15.tsv`.

**Proposed:**
- Add the 245 models flagged `concordant`: 51 cell lines (mostly Maris
  neuroblastoma) and 194 PDXs.
  - `cell_line_id`: the OpenPedCan `sample_id`, prefixed by cohort.
  - `source`: `OpenPedCan v15 (<cohort>)`.
- Keep `intermediate` and `discordant` models out of the default views, or show
  them flagged.

| field | value |
|---|---|
| `tumor_concordance_score` | `concordance` (Spearman vs own-group tumour centroid) |
| `tumor_concordance_flag` | `flag` |

PDXs are not cell lines. If `cell_line` should stay cell-line-only, add a
`model_type` column or a separate table rather than mixing them silently.

OncoTree codes come from `configs/ccdi/diagnosis_map.tsv`. They were curated
by hand and not checked live against oncotree.info.

## 2. Existing NPxP cell lines

No change. None of the 16 has expression in OpenPedCan, and none is a
pediatric model. `tumor_concordance_score` stays null, and the flag stays
`unassessed` rather than becoming `no_cohort`, because no comparison was
attempted. Scoring them requires DepMap expression, which the backend
workspace has.

## 3. Mechanism hypotheses

Already delivered as a pipeline-produced page (same precedent as the original
export):
- `api/mechanism_hypotheses.json`, with `pediatric_context` per row and a
  top-level gene table;
- `js/mechanism.js`, updated to show it.

If the backend regenerates this file, it must call `build-ccdi --step
annotate` after building the hypotheses, or the pediatric blocks will
disappear.

## 4. Opportunity Ranking burden weighting

Not done: SEER was unreachable. When it is, add the burden columns
(incidence, 5-year survival, ages 0–19) alongside the existing score. Do not
change the score formula without a separately documented methodology change.
