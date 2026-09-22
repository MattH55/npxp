# NPxP — Non-Pharm × Pharm

An open, evidence-graded database modeling interactions between non-pharmacological interventions (dietary/metabolic, thermal, hypoxic, and similar modifiers) and pharmacological agents, at the cellular level, across cancer cell lines.

Part of the [Open Source Medicine Foundation](https://opensourcemed.info) portfolio. Hosted at [npxp.opensourcemed.info](https://npxp.opensourcemed.info).

See `BUILD_SPEC.md` for the full data model, ingestion sources, API surface, and page structure.

## Repository layout (GitHub Pages deployment)

- `*.html`, `js/`, `css/` — the static frontend (vanilla JS + Tailwind CDN).
- `api/` — static JSON snapshot of the REST API responses; the frontend reads
  these files because GitHub Pages cannot run the FastAPI backend.
- `js/config.js` — sets `STATIC_DATA = true`; locally the same frontend runs
  against the live API instead (`STATIC_DATA = false`).
- `js/scoring-client.js` — browser port of the Step 1 synergy engine
  (Bliss/HSA/Loewe/ZIP), used by the scoring page in the static build.
- `CNAME` — custom domain `npxp.opensourcemed.info`.

The backend, ingestion pipeline, curation tooling, and tests live in the
development workspace (`synlethality/` package). After any data or frontend
change, regenerate this snapshot with:

```bash
python scripts/export_static.py --repo <path-to-this-repo>
```

No quantitative field is ever fabricated: `stress_signature_score.score` and
`interaction_effect.combined_effect_metric` are null until computed from raw
data by the ingestion pipeline.
