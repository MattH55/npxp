"""Build the real, mechanism-cited modifier x drug hypothesis table.

For each modifier in configs/npi_mechanisms.yaml:
  - induced_hr_deficiency: look up real synthetic-lethal partners of the
    modifier's lesion gene(s) in SynLethDB, then real approved anti-
    neoplastic drugs targeting each partner in DGIdb.
  - pathway_suppression: real approved anti-neoplastic drugs directly
    targeting the literature-recommended gene(s), from DGIdb.
  - immune_mobilization: real approved immunotherapy (checkpoint inhibitor)
    drugs, from DGIdb.

Every row in the output traces to: the modifier's own PMID (why this
mechanism is real), plus, for SL-derived rows, the SL pair's own PMID (why
that gene pair is a real synthetic-lethal partner) -- two independent real
citations, never a numeric score standing in for either.

    python scripts/build_mechanism_hypotheses.py --out out/mechanism_hypotheses.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone

import yaml

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

from npi_pharma.mechanism.dgidb import (
    load_gene_to_drugs,
    real_anticancer_drugs,
    real_immunotherapy_drugs,
)
from npi_pharma.mechanism.synlethdb import load_sl_pairs, sl_partners

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_PATH = os.path.join(ROOT, "configs", "npi_mechanisms.yaml")

# Canonical immune-checkpoint genes -- see dgidb.py's DIRECT_ACTION_TYPES
# docstring for why a plain "immunotherapy" flag alone is too broad
# (catches immunosuppressants like sirolimus).
CHECKPOINT_GENES = ["PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT"]

#: Real SL evidence-type quality tiers (highest first), per SynLethDB's own
#: rel_source vocabulary -- a direct CRISPR/CRISPRi screen is a measured
#: double-knockout result; Text Mining is an NLP-extracted literature
#: claim, the weakest real evidence type SynLethDB carries (not a
#: judgment call invented here -- SynLethDB's own methods paper describes
#: the same hierarchy of experimental vs. literature-derived evidence).
SL_EVIDENCE_TIER = {
    "CRISPR/CRISPRi": 3, "High Throughput": 2, "Low Throughput": 2,
    "Synlethality": 2, "Decipher": 2, "Computational Prediction": 1,
    "Text Mining": 0,
}


def _evidence_tier(evidence_type: str | None) -> int:
    if not evidence_type:
        return 0
    # rel_source can combine multiple types with ';' -- take the best one.
    return max((SL_EVIDENCE_TIER.get(t.strip(), 0) for t in evidence_type.split(";")), default=0)


def build_hr_deficiency_rows(modifier_id: str, spec: dict, sl_index, gene_index) -> list[dict]:
    rows = []
    for lesion_gene in spec["lesion_genes"]:
        for partner in sl_partners(lesion_gene, sl_index):
            drugs = real_anticancer_drugs(partner["partner_gene"], gene_index)
            if not drugs:
                continue
            for drug in drugs:
                rows.append({
                    "modifier_id": modifier_id,
                    "modifier_label": spec["label"],
                    "category": "induced_hr_deficiency",
                    "lesion_gene": lesion_gene,
                    "sl_partner_gene": partner["partner_gene"],
                    "sl_evidence_type": partner["evidence_type"],
                    "sl_evidence_tier": _evidence_tier(partner["evidence_type"]),
                    "sl_pubmed_id": partner["pubmed_id"],
                    "drug_name": drug["drug_name"],
                    "drug_target_gene": partner["partner_gene"],
                    "drug_n_sources": drug["n_sources"],
                    "modifier_citation": spec["citation"],
                    "supporting_citations": spec.get("supporting_citations", []),
                    "mechanism": spec["mechanism"].strip(),
                })
    return rows


def build_pathway_suppression_rows(modifier_id: str, spec: dict, gene_index) -> list[dict]:
    rows = []
    for gene in spec["recommended_target_genes"]:
        for drug in real_anticancer_drugs(gene, gene_index):
            rows.append({
                "modifier_id": modifier_id,
                "modifier_label": spec["label"],
                "category": "pathway_suppression",
                "lesion_gene": gene,
                "sl_partner_gene": None,
                "sl_evidence_type": None,
                "sl_evidence_tier": None,
                "sl_pubmed_id": None,
                "drug_name": drug["drug_name"],
                "drug_target_gene": gene,
                "drug_n_sources": drug["n_sources"],
                "modifier_citation": spec["citation"],
                "supporting_citations": spec.get("supporting_citations", []),
                "mechanism": spec["mechanism"].strip(),
            })
    return rows


def build_immune_mobilization_rows(modifier_id: str, spec: dict, gene_index) -> list[dict]:
    rows = []
    for gene in CHECKPOINT_GENES:
        for drug in real_anticancer_drugs(gene, gene_index):
            rows.append({
                "modifier_id": modifier_id,
                "modifier_label": spec["label"],
                "category": "immune_mobilization",
                "lesion_gene": gene,
                "sl_partner_gene": None,
                "sl_evidence_type": None,
                "sl_evidence_tier": None,
                "sl_pubmed_id": None,
                "drug_name": drug["drug_name"],
                "drug_target_gene": gene,
                "drug_n_sources": drug["n_sources"],
                "modifier_citation": spec["citation"],
                "supporting_citations": spec.get("supporting_citations", []),
                "mechanism": spec["mechanism"].strip(),
            })
    return rows


def _attach_clinical_trials(all_rows: list[dict], trials_path: str) -> None:
    """Join the real ClinicalTrials.gov cross-check (scripts/fetch_clinical_trials.py)
    onto matching rows, keyed exactly as that script queried them
    (f"{modifier_id}::{drug_name}") -- so a row with no entry in the file
    (not run for this modifier, or genuinely zero real trials found) is
    left with an empty list, never fabricated as a gap."""
    if not os.path.exists(trials_path):
        print(f"  (no clinical trials file at {trials_path} -- skipping join)")
        return
    with open(trials_path, encoding="utf-8") as fh:
        ct = json.load(fh)
    for r in all_rows:
        key = f"{r['modifier_id']}::{r['drug_name']}"
        entry = ct["results"].get(key)
        r["clinical_trials"] = entry["trials"] if entry else []
        r["has_clinical_trial"] = bool(r["clinical_trials"])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=os.path.join(ROOT, "out", "mechanism_hypotheses.json"))
    parser.add_argument("--clinical-trials", default=os.path.join(ROOT, "out", "clinical_trials.json"))
    args = parser.parse_args()

    with open(CONFIG_PATH, encoding="utf-8") as fh:
        config = yaml.safe_load(fh)

    print("Loading real SynLethDB SL pairs and DGIdb drug-gene interactions")
    sl_index = load_sl_pairs()
    gene_index = load_gene_to_drugs()
    print(f"  {len(sl_index)} genes in SynLethDB, {len(gene_index)} genes in DGIdb")

    all_rows = []
    for modifier_id, spec in config["modifiers"].items():
        category = spec["category"]
        print(f"Resolving {modifier_id} ({category})")
        if category == "induced_hr_deficiency":
            rows = build_hr_deficiency_rows(modifier_id, spec, sl_index, gene_index)
        elif category == "pathway_suppression":
            rows = build_pathway_suppression_rows(modifier_id, spec, gene_index)
        elif category == "immune_mobilization":
            rows = build_immune_mobilization_rows(modifier_id, spec, gene_index)
        else:
            raise ValueError(f"{modifier_id}: unknown category {category!r}")
        print(f"  {len(rows)} real, cited hypothesis rows")
        all_rows.extend(rows)

    # Rank strongest-evidence-first: real SL evidence tier (direct CRISPR
    # screen > literature-curated > computational > text-mined), then how
    # many independent DGIdb sources corroborate the drug-target claim.
    # pathway_suppression/immune_mobilization rows have no SL tier (the
    # mechanism is literature-direct, not SL-derived) -- treated as tier 2
    # so they rank alongside literature-curated SL evidence, not above or
    # below it by construction.
    def _direct_support(row: dict) -> list[dict]:
        """Real supporting citations that specifically tested THIS row's
        drug-target gene, not just the same modifier generically -- e.g.
        the heat+PARP-inhibitor follow-up paper covers PARP1, not every
        other SL partner of BRCA2 a drug might be found for."""
        return [
            c for c in row["supporting_citations"]
            if row["drug_target_gene"] in c.get("tested_target_genes", [])
        ]

    for r in all_rows:
        r["direct_experimental_support"] = _direct_support(r)
        r["has_direct_experimental_support"] = bool(r["direct_experimental_support"])

    all_rows.sort(
        key=lambda r: (r["has_direct_experimental_support"],
                       r["sl_evidence_tier"] if r["sl_evidence_tier"] is not None else 2,
                       r["drug_n_sources"]),
        reverse=True,
    )

    print("Joining real ClinicalTrials.gov cross-check (scripts/fetch_clinical_trials.py)")
    _attach_clinical_trials(all_rows, args.clinical_trials)
    n_with_trial = sum(1 for r in all_rows if r.get("has_clinical_trial"))
    print(f"  {n_with_trial}/{len(all_rows)} rows have a matching real registered trial")

    out = {
        "built": datetime.now(timezone.utc).isoformat(),
        "method": (
            "Mechanism-cited hypothesis generation, NOT a calibrated synergy score. "
            "Each row traces to: (1) a real PMID for why the modifier induces the "
            "stated gene/pathway state, (2) for induced_hr_deficiency rows, a real "
            "SynLethDB-curated synthetic-lethal pair (with its own PMID where "
            "available) linking the induced lesion gene to the drug's target gene, "
            "and (3) a real DGIdb drug-gene interaction claim, counted by how many "
            "independent source databases report it. This project tested whether "
            "expression-signature similarity/composition predicts real measured "
            "synergy (both a pathway-feature model here and a much larger "
            "740k-pair DrugComb test) and found no signal either time -- see "
            "pipeline/docs/validation_drugcomb.md. This table does not attempt "
            "that; it is real literature + real curated databases, not a model. "
            "Each row also carries `clinical_trials`: a real ClinicalTrials.gov "
            "cross-check (scripts/fetch_clinical_trials.py) for whether this "
            "exact modifier+drug combination has already been registered as a "
            "trial -- an empty list means no matching trial was found (a "
            "candidate white-space gap), not that none exists; the search is "
            "not exhaustive (see that script's own limitations note)."
        ),
        "n_rows": len(all_rows),
        "rows": all_rows,
    }
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nWrote {len(all_rows)} total rows to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
