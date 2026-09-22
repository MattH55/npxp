# NPxP Prediction Methodology

How NPxP predicts drug × non-pharmacological modifier interaction effects for combinations that have never been directly tested — in three stages, ordered by what's actually available today versus what has to be earned with accumulating data. Each stage produces output at a specific `evidence_tier`, so a prediction is never presented with more confidence than its method actually supports.

This doc is the technical backbone of the Prediction Layer and Mechanistic Bridging module described in `BUILD_SPEC.md` — read that first for the surrounding schema (`stress_signature_score`, `drug_resistance_mechanism`, `interaction_effect`).

---

## Stage 1 — Signature-overlap scoring (available now, no model required)

**When it applies:** any drug × modifier × cell-line triple where both agents have an induced gene-expression signature — i.e. almost everything in the library from day one, since drugs come from LINCS L1000 and modifiers come from the GEO studies already ingested into `stress_signature_score`.

**Method:** compare the two signatures — significantly changed genes, enriched pathways/gene sets (MSigDB Hallmark, FerrDb, etc.) — and score their overlap and directional concordance. This is the "co-gene/GS" approach, benchmarked against the DREAM consortium's gold-standard synergy dataset and shown to outperform dedicated DREAM-challenge submissions (PC-index 0.663), despite being computationally simple. The purpose-built `ccmap` (Combination Connectivity Mapping) R package implements the combination-level version of this directly — correlating a combined signature against a query — and is the recommended starting implementation rather than writing the overlap logic from scratch.

**Directionality matters as much as overlap:** a modifier that overlaps a drug's resistance pathway and *suppresses* it predicts synergy; overlap with the pathway being *reinforced* instead predicts antagonism (the pattern already seen with HSP70/hyperthermia and doxorubicin resistance). Don't score overlap alone — always resolve direction against the `drug_resistance_mechanism` entry before assigning `interaction_type`.

**Output:** `interaction_effect` rows with `evidence_tier = tier_3_mechanism_only`. This is exactly the Mechanistic Bridging module's output — no death/viability readout has been measured for the combination itself, only signature-level plausibility. Never let this tier read as more confident than it is in the UI.

**Why this comes first:** it requires no training data, no labeled ground truth, and no model to build or maintain — it's pure computation on signatures you already have. This is what makes Coverage Gap analysis and the Mechanistic Bridging candidate list functional from launch, before any dose-response data has accumulated.

---

## Stage 2 — Classical dose-response synergy quantification (as Tier 1 data accumulates)

**When it applies:** once an actual dose-response matrix exists for a specific drug × modifier × cell-line combination — e.g. BHB-concentration × drug-dose, or hyperthermia-temperature × drug-dose, measured in the same study (the 7-cell-line BHB/chemo-radiosensitivity study, the PEMF/temozolomide study, etc.).

**Method:** apply Bliss independence, Loewe additivity, HSA, and/or ZIP scoring — the same standard models SynergyFinder and DrugComb use for drug-drug pairs. None of these require the second factor to be a drug; they operate on any dose-response surface, so they apply unchanged to modifier × drug matrices.

**Output:** `interaction_effect` rows with `evidence_tier = tier_1_direct`, carrying an actual `combined_effect_metric`. These rows are the ground truth — both for direct display to users and as the labeled training data Stage 3 depends on.

**Practical note:** this stage is gated by data availability, not by build complexity — the math is simple and should be implemented early, but its output volume grows only as fast as Tier 1 studies are ingested (see `BUILD_SPEC.md`'s Modifier Scope table and Coverage Gap analysis for what's actually available today).

---

## Stage 3 — MARSY fine-tuning (once there's enough Tier 1 ground truth)

**When it applies:** after Stage 2 has produced a meaningful volume of labeled Tier 1 drug × modifier combinations across enough cell lines and modifier types to fine-tune a model responsibly — not before. There's no fixed threshold to commit to in advance; treat this as a go/no-go decision revisited periodically as Tier 1 data grows, not a fixed milestone date.

**Base model:** [MARSY](https://github.com/Emad-COMBINE-lab/MARSY) (Multitask drug pAiR SynergY) — publicly released, with code and cleaned training data on GitHub, and pre-published predictions for ~133,000 drug-pair × cell-line combinations as a community resource. Architecturally, MARSY represents each agent via **(a)** the cell line's baseline gene expression and **(b)** the agent's induced differential-expression signature — not chemical structure — which is exactly the representation your modifiers already have via `stress_signature_score`. That compatibility is what makes MARSY the right base model rather than a chemical-structure-based alternative like DeepSynergy.

**The honest limitation:** MARSY's published weights were trained and validated exclusively on drug-drug pairs. Applying them unmodified to a drug × modifier pair is an out-of-distribution use with no validated accuracy guarantee, even though the input format lines up. **This is transfer learning, not a plug-and-play oracle:** use the public weights as an initialization, then fine-tune on NPxP's own accumulated Tier 1 drug × modifier ground truth from Stage 2, rather than either (a) using the weights unmodified, or (b) training an equivalent architecture from random initialization. Fine-tuning from a pretrained starting point needs meaningfully less labeled data than training from scratch — but it still needs real Tier 1 data to fine-tune against; a thin dataset produces a thin fine-tune, not a shortcut around the data requirement.

**Output:** `interaction_effect` rows with `evidence_tier = tier_2b_model_predicted`, each one storing the model version/fine-tuning run ID it came from, so predictions can be invalidated and regenerated as the training data or model improves rather than silently overwritten.

**Until this stage is reached:** the Prediction Layer and Mechanistic Bridging outputs run entirely on Stages 1 and 2. That's a legitimate MVP state, not a placeholder — it should ship fully functional on its own, with Stage 3 added later as a genuine upgrade rather than something the launch depends on.

---

## Evidence tier summary

| Stage | Method | Evidence tier | Data/model prerequisite |
|---|---|---|---|
| 1 | Signature-overlap scoring (co-gene/GS, `ccmap`) | `tier_3_mechanism_only` | None — signatures already in hand |
| 2 | Bliss/Loewe/HSA/ZIP on real dose-response data | `tier_1_direct` | A study testing the combination directly |
| 3 | MARSY fine-tuned on NPxP's own Tier 1 data | `tier_2b_model_predicted` | Sufficient accumulated Tier 1 volume from Stage 2 |

Separately, `tier_2_inferred` (manual join of separate single-factor studies via shared cell line and plausible mechanism, without either a model or a formal signature-overlap computation) remains a distinct, human-curated category — not produced by any of the three stages above.
