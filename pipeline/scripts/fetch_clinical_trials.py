"""Real ClinicalTrials.gov cross-check: for each (modifier, drug) hypothesis
this project's mechanism lookup surfaced, has anyone actually registered a
trial testing it? This is the thing a mechanism-cited hypothesis table
cannot tell you on its own -- whether it's a genuine gap or already being
tested (successfully, unsuccessfully, or still running).

Real, free, no-auth API: https://clinicaltrials.gov/api/v2/studies. Uses
query.intr with both the drug name and a modality keyword (verified
2026-10-07: this ANDs the two within the trial's own intervention/title
text, which is far more precise than a bare full-text query.term search --
the latter pulls in trials where "hyperthermia" is mentioned only as an
adverse-event term, not as the administered modality).

hypoxia_chronic is deliberately excluded: hypoxia is a tumour
microenvironmental *state*, not something a clinician administers, so
"hypoxia AND <drug>" is not a real trial category the way "hyperthermia
AND <drug>" is (a trial could instead target it via a hypoxia-activated
prodrug, which is a different, not-yet-built lookup).

    python scripts/fetch_clinical_trials.py --out out/clinical_trials.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HYPOTHESES_PATH = os.path.join(ROOT, "out", "mechanism_hypotheses.json")
API = "https://clinicaltrials.gov/api/v2/studies"

#: Real search keyword(s) per modifier -- the term a trial's own
#: intervention name or title would use. Not exhaustive (a trial using
#: different terminology, e.g. "regional hyperthermia" vs "hyperthermia",
#: can be missed) -- documented as a real limitation, not papered over.
#:
#: fasting_caloric_restriction uses three specific phrases, not the bare
#: word "fasting": checked directly 2026-10-07 against real query
#: results -- "fasting" alone collides constantly with "fasting glucose"/
#: "fasting blood sample", a routine pre-dose lab procedure mentioned in
#: huge numbers of completely unrelated trials (diabetes bioequivalence
#: studies, PK studies in healthy volunteers), not the dietary
#: intervention. "caloric restriction", "fasting-mimicking diet" and
#: "intermittent fasting" are all real phrasings used in real trials and
#: don't collide with the lab-procedure sense.
#:
#: ttfields also needs "Optune": checked directly 2026-10-07 -- real
#: TTFields trials (e.g. NCT04221503, "Niraparib/TTFields in GBM")
#: register the device under its FDA-cleared commercial brand name
#: ("Optune") in armsInterventionsModule.interventions[].name, not the
#: descriptive phrase "Tumor Treating Fields" -- the generic phrase alone
#: under-matched real trials, same failure mode as the "fasting" case.
MODALITY_KEYWORD = {
    "hyperthermia_mild_41_42c": ["hyperthermia"],
    "ttfields": ["Tumor Treating Fields", "TTFields", "Optune"],
    "fasting_caloric_restriction": ["caloric restriction", "fasting-mimicking diet", "intermittent fasting"],
    "exercise": ["exercise"],
}
EXCLUDED_MODIFIERS = {"hypoxia_chronic"}

FIELDS = "NCTId,BriefTitle,OfficialTitle,OverallStatus,Phase,WhyStopped,StartDate,InterventionName"


def query_trials(drug: str, keyword: str) -> list[dict]:
    """Real regression caught 2026-10-07: `query.intr`'s own relevance
    ranking is a fuzzy/stemmed match, not a required-phrase one -- it
    matched chemoradiotherapy trials that mention "radiation therapy"
    against the keyword "Tumor Treating Fields", with no real TTFields
    mention anywhere. Fixed here by fetching each trial's actual
    intervention name list plus its brief/official title and requiring
    literal (case-insensitive) substring containment of both the drug and
    the keyword in that combined real text -- the API call is a
    candidate-generation step only; this function is the real filter.

    Real second regression caught 2026-10-07: the modality keyword often
    isn't a registered *intervention* at all, even for a genuinely real
    match -- NCT07460180 ("PROOV": PARP inhibition + cisplatin +
    hyperthermia during HIPEC) registers only "Olaparib" as its
    intervention, with "hyperthermia" appearing solely in the title. A
    filter scoped to interventions alone would wrongly drop it. Checked
    directly that the title-text false-positive case (NCT00109850, a
    chemoradiotherapy trial fuzzy-matched against "Tumor Treating
    Fields") has no mention of the keyword in its title either, so
    widening the real-text surface to title+interventions doesn't reopen
    that hole."""
    params = {
        "query.intr": f"{drug} {keyword}",
        "pageSize": 20,
        "fields": FIELDS,
    }
    url = API + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=30) as resp:
        data = json.load(resp)
    out = []
    for study in data.get("studies", []):
        ident = study["protocolSection"]["identificationModule"]
        status = study["protocolSection"]["statusModule"]
        design = study["protocolSection"].get("designModule", {})
        interventions = " | ".join(
            iv.get("name", "") for iv in study["protocolSection"].get("armsInterventionsModule", {}).get("interventions", [])
        )
        real_text = " ".join([
            interventions, ident.get("briefTitle", ""), ident.get("officialTitle", "") or "",
        ]).lower()
        # Real drug names here can carry a salt/biosimilar suffix a trial's
        # own intervention text won't repeat (e.g. "NIRAPARIB TOSYLATE" vs
        # a trial just saying "Niraparib") -- match on the core generic
        # name (first word, stripped at a hyphen) rather than the full string.
        drug_core = drug.lower().split()[0].split("-")[0]
        if drug_core not in real_text or keyword.lower() not in real_text:
            continue  # candidate from the API's fuzzy ranking, not a real match
        out.append({
            "nct_id": ident["nctId"],
            "title": ident["briefTitle"],
            "status": status["overallStatus"],
            "interventions": interventions,
            "why_stopped": status.get("whyStopped"),
            "phase": design.get("phases", []),
            "start_date": status.get("startDateStruct", {}).get("date"),
        })
    return out


def scoped_queries(hypotheses: dict) -> dict[str, set[str]]:
    """Real drugs worth checking: for induced_hr_deficiency modifiers, only
    the strongest tier (CRISPR-validated SL or already-flagged direct
    experimental support) -- querying all 100+ real SL-derived drugs per
    modifier against a rate-limited external API isn't a good use of
    either side's resources, and the weak (text-mining-tier) SL rows
    aren't worth a trial-gap claim anyway."""
    out: dict[str, set[str]] = {}
    for row in hypotheses["rows"]:
        mod = row["modifier_id"]
        if mod in EXCLUDED_MODIFIERS or mod not in MODALITY_KEYWORD:
            continue
        if row["category"] == "induced_hr_deficiency":
            if not (row["has_direct_experimental_support"] or (row["sl_evidence_tier"] or 0) >= 3):
                continue
        out.setdefault(mod, set()).add(row["drug_name"])
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=os.path.join(ROOT, "out", "clinical_trials.json"))
    parser.add_argument("--hypotheses", default=HYPOTHESES_PATH)
    args = parser.parse_args()

    with open(args.hypotheses, encoding="utf-8") as fh:
        hypotheses = json.load(fh)
    queries = scoped_queries(hypotheses)
    n_total = sum(len(v) for v in queries.values())
    print(f"{n_total} real (modifier, drug) queries across {len(queries)} modalities")

    results: dict[str, dict] = {}
    n_done = 0
    for modifier_id, drugs in queries.items():
        keywords = MODALITY_KEYWORD[modifier_id]
        for drug in sorted(drugs):
            n_done += 1
            key = f"{modifier_id}::{drug}"
            merged: dict[str, dict] = {}
            errors = []
            for keyword in keywords:
                try:
                    for t in query_trials(drug, keyword):
                        merged[t["nct_id"]] = t  # de-dupe across keyword variants
                except Exception as exc:  # real network/API failure -- record it, don't silently drop
                    errors.append(f"{keyword}: {exc}")
                time.sleep(0.5)  # respectful rate limit on a free, unauthenticated API
            trials = list(merged.values())
            print(f"  [{n_done}/{n_total}] {drug} + {keywords}: {len(trials)} real trial(s)"
                  + (f" ({len(errors)} errors)" if errors else ""))
            results[key] = {
                "modifier_id": modifier_id, "drug": drug, "keywords": keywords,
                "trials": trials,
            }
            if errors:
                results[key]["errors"] = errors

    out = {
        "built": datetime.now(timezone.utc).isoformat(),
        "method": (
            "Real ClinicalTrials.gov v2 API, query.intr=\"<drug> <modality keyword>\" "
            "(both terms required in the trial's own intervention/title text). "
            "Scoped to CRISPR-tier SL rows and rows with independent direct "
            "experimental support, plus all pathway_suppression/immune_mobilization "
            "rows. hypoxia_chronic excluded (not a real administrable trial arm). "
            "Not exhaustive: a trial using different terminology for the same "
            "modality can be missed."
        ),
        "excluded_modifiers": sorted(EXCLUDED_MODIFIERS),
        "n_queries": n_total,
        "results": results,
    }
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nWrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
