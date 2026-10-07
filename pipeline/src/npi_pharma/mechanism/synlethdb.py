"""Real synthetic-lethality lookups against SynLethDB 3.0's curated human
SL table (data/raw/synlethdb/Human.SL.detailed.tsv, fetched by
scripts/fetch_mechanism_data.py). Every row is a real, curated
(non-computational) SL pair with a real PubMed citation and evidence type.
"""
from __future__ import annotations

import csv
import os
from collections import defaultdict
from dataclasses import dataclass

DEFAULT_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))),
    "data", "raw", "synlethdb", "Human.SL.detailed.tsv",
)


@dataclass
class SLPair:
    gene_a: str
    gene_b: str
    evidence_type: str
    pubmed_id: str
    cell_line: str


def load_sl_pairs(path: str = DEFAULT_PATH) -> dict[str, list[SLPair]]:
    """{gene_symbol -> [SLPair, ...]}, indexed on both sides of every real
    pair (SL is symmetric; the source file lists each pair once)."""
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"{path} not found. Run `python scripts/fetch_mechanism_data.py` first."
        )
    by_gene: dict[str, list[SLPair]] = defaultdict(list)
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get(":TYPE") != "Gene_SL_Gene":
                continue
            a, b = row["x_name"].strip(), row["y_name"].strip()
            if not a or not b:
                continue
            pair = SLPair(
                gene_a=a, gene_b=b,
                evidence_type=row.get("rel_source", ""),
                pubmed_id=row.get("pubmed_id", ""),
                cell_line=row.get("cell_line", ""),
            )
            by_gene[a].append(pair)
            by_gene[b].append(pair)
    return dict(by_gene)


def sl_partners(gene: str, sl_index: dict[str, list[SLPair]]) -> list[dict]:
    """Real SL partner genes of `gene`, each with its own real citation.
    One row per distinct (partner, pubmed_id) -- a partner reported by
    multiple independent studies keeps each citation, not collapsed."""
    out = []
    seen = set()
    for pair in sl_index.get(gene, []):
        partner = pair.gene_b if pair.gene_a == gene else pair.gene_a
        key = (partner, pair.pubmed_id)
        if key in seen:
            continue
        seen.add(key)
        out.append({
            "partner_gene": partner,
            "evidence_type": pair.evidence_type,
            "pubmed_id": pair.pubmed_id,
            "cell_line": pair.cell_line,
        })
    return out
