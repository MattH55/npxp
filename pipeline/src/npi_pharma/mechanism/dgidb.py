"""Real drug-gene interaction lookups against DGIdb 5.0's bulk interactions
export (data/raw/dgidb/interactions.tsv, fetched by
scripts/fetch_mechanism_data.py).
"""
from __future__ import annotations

import csv
import os
from collections import defaultdict
from dataclasses import dataclass

DEFAULT_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))),
    "data", "raw", "dgidb", "interactions.tsv",
)


@dataclass
class DrugHit:
    drug_name: str
    gene: str
    interaction_type: str
    source_db: str
    approved: bool
    immunotherapy: bool
    anti_neoplastic: bool


def load_gene_to_drugs(path: str = DEFAULT_PATH) -> dict[str, list[DrugHit]]:
    """{gene_symbol -> [DrugHit, ...]}, every real DGIdb claim (deduplication
    is the caller's job -- a gene commonly has one real drug hit reported by
    several source databases, which is corroboration, not noise)."""
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"{path} not found. Run `python scripts/fetch_mechanism_data.py` first."
        )
    out: dict[str, list[DrugHit]] = defaultdict(list)
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene = (row.get("gene_name") or "").strip()
            drug = (row.get("drug_name") or "").strip()
            if not gene or not drug:
                continue
            out[gene].append(DrugHit(
                drug_name=drug, gene=gene,
                interaction_type=row.get("interaction_type", "") or "",
                source_db=row.get("interaction_source_db_name", ""),
                approved=row.get("approved") == "TRUE",
                immunotherapy=row.get("immunotherapy") == "TRUE",
                anti_neoplastic=row.get("anti_neoplastic") == "TRUE",
            ))
    return dict(out)


#: Real "this drug does something to this gene" verbs, as opposed to a bare
#: biomarker/association claim (DGIdb's interaction_type is NULL for a large
#: share of rows that only report co-occurrence in a biomarker study -- e.g.
#: CTLA4 "associated with" cyclosporine or prednisone, which are not CTLA4-
#: targeting drugs). Checked directly 2026-10-07: without this filter,
#: PDCD1/CD274/CTLA4 lookups return mostly unrelated chemotherapies and
#: immunosuppressants; with it, exactly the real checkpoint inhibitors.
DIRECT_ACTION_TYPES = {
    "inhibitor", "antagonist", "blocker", "agonist", "activator",
    "inducer", "modulator", "partial agonist", "allosteric modulator",
}


def real_anticancer_drugs(
    gene: str, gene_index: dict[str, list[DrugHit]], require_approved: bool = True,
    require_direct_action: bool = True,
) -> list[dict]:
    """Real, distinct approved anti-neoplastic drugs targeting `gene`, each
    with the set of source databases that independently reported it.
    `require_direct_action`: keep a drug only if at least one real source
    reports a direct-action interaction_type for it, not only a bare
    biomarker-association claim (see DIRECT_ACTION_TYPES)."""
    by_drug: dict[str, dict] = {}
    for hit in gene_index.get(gene, []):
        if not hit.anti_neoplastic:
            continue
        if require_approved and not hit.approved:
            continue
        entry = by_drug.setdefault(hit.drug_name, {
            "drug_name": hit.drug_name, "gene": gene, "approved": hit.approved,
            "interaction_types": set(), "source_dbs": set(), "has_direct_action": False,
        })
        if hit.interaction_type:
            entry["interaction_types"].add(hit.interaction_type)
            if hit.interaction_type.lower() in DIRECT_ACTION_TYPES:
                entry["has_direct_action"] = True
        entry["source_dbs"].add(hit.source_db)
    out = []
    for entry in by_drug.values():
        if require_direct_action and not entry["has_direct_action"]:
            continue
        entry["interaction_types"] = sorted(entry["interaction_types"])
        entry["source_dbs"] = sorted(entry["source_dbs"])
        entry["n_sources"] = len(entry["source_dbs"])
        out.append(entry)
    out.sort(key=lambda e: e["n_sources"], reverse=True)
    return out


def real_immunotherapy_drugs(gene_index: dict[str, list[DrugHit]]) -> list[dict]:
    """Every real approved immunotherapy drug in DGIdb, regardless of gene
    (used for the exercise/immune_mobilization pairing, which isn't a
    single-gene lookup)."""
    by_drug: dict[str, dict] = {}
    for hits in gene_index.values():
        for hit in hits:
            if not (hit.immunotherapy and hit.approved):
                continue
            entry = by_drug.setdefault(hit.drug_name, {
                "drug_name": hit.drug_name, "genes": set(), "source_dbs": set(),
            })
            entry["genes"].add(hit.gene)
            entry["source_dbs"].add(hit.source_db)
    out = []
    for entry in by_drug.values():
        entry["genes"] = sorted(entry["genes"])
        entry["source_dbs"] = sorted(entry["source_dbs"])
        entry["n_sources"] = len(entry["source_dbs"])
        out.append(entry)
    out.sort(key=lambda e: e["n_sources"], reverse=True)
    return out
