"""Unit tests for npi_pharma.mechanism against small synthetic fixtures --
never the real downloaded SynLethDB/DGIdb files, same convention as the
rest of this project's tests (see test_ingest.py)."""
import csv
import importlib.util
import json
import os
import sys

import pytest

from npi_pharma.mechanism.dgidb import DIRECT_ACTION_TYPES, load_gene_to_drugs, real_anticancer_drugs, real_immunotherapy_drugs
from npi_pharma.mechanism.synlethdb import load_sl_pairs, sl_partners

_SCRIPT_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts", "build_mechanism_hypotheses.py"
)
_spec = importlib.util.spec_from_file_location("build_mechanism_hypotheses", _SCRIPT_PATH)
build_mechanism_hypotheses = importlib.util.module_from_spec(_spec)
sys.modules["build_mechanism_hypotheses"] = build_mechanism_hypotheses
_spec.loader.exec_module(build_mechanism_hypotheses)
_attach_clinical_trials = build_mechanism_hypotheses._attach_clinical_trials


def _write_sl_tsv(path, rows):
    header = ["x:START_ID", "x_type", "x_name", "x_source", "y:END_ID", "y_type", "y_name",
              "y_source", "relation", ":TYPE", "rel_source", "edge_index", "cell_line",
              "pubmed_id", "cancer"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for a, b, evidence, pmid in rows:
            w.writerow([1, "Gene", a, "NCBI", 2, "Gene", b, "NCBI", "SL", "Gene_SL_Gene",
                        evidence, 999, "TBD", pmid, ""])


def _write_dgidb_tsv(path, rows):
    header = ["gene_claim_name", "gene_concept_id", "gene_name", "interaction_source_db_name",
              "interaction_source_db_version", "interaction_type", "interaction_score",
              "drug_claim_name", "drug_concept_id", "drug_name", "approved", "immunotherapy",
              "anti_neoplastic"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for gene, source, itype, drug, approved, immuno, antineo in rows:
            w.writerow([gene, "hgnc:0", gene, source, "1", itype, "0.5", drug, "rxcui:0",
                        drug, str(approved).upper(), str(immuno).upper(), str(antineo).upper()])


@pytest.fixture
def sl_tsv(tmp_path):
    path = tmp_path / "sl.tsv"
    _write_sl_tsv(path, [
        ("BRCA2", "PARP1", "CRISPR/CRISPRi", "21555554"),
        ("BRCA2", "PARP2", "Low Throughput", "28628639"),
        ("ATM", "ATR", "Text Mining", ""),
    ])
    return path


@pytest.fixture
def dgidb_tsv(tmp_path):
    path = tmp_path / "interactions.tsv"
    _write_dgidb_tsv(path, [
        ("PARP1", "ChEMBL", "inhibitor", "OLAPARIB", True, False, True),
        ("PARP1", "CIViC", "NULL", "OLAPARIB", True, False, True),  # corroborating, no type
        ("PARP1", "TTD", "inhibitor", "NIRAPARIB", True, False, True),
        ("PARP1", "DTC", "NULL", "EXPERIMENTAL-X", False, False, True),  # not approved
        ("CTLA4", "CIViC", "NULL", "PREDNISONE", True, False, False),  # biomarker noise, not anti-neoplastic
        ("CTLA4", "ChEMBL", "inhibitor", "IPILIMUMAB", True, True, True),
    ])
    return path


class TestSynLethDB:
    def test_load_and_index_both_directions(self, sl_tsv):
        idx = load_sl_pairs(str(sl_tsv))
        assert "BRCA2" in idx and "PARP1" in idx  # indexed on both sides

    def test_sl_partners_real_pair(self, sl_tsv):
        idx = load_sl_pairs(str(sl_tsv))
        partners = sl_partners("BRCA2", idx)
        names = {p["partner_gene"] for p in partners}
        assert names == {"PARP1", "PARP2"}

    def test_sl_partners_carries_citation(self, sl_tsv):
        idx = load_sl_pairs(str(sl_tsv))
        partners = sl_partners("BRCA2", idx)
        parp1 = next(p for p in partners if p["partner_gene"] == "PARP1")
        assert parp1["pubmed_id"] == "21555554"
        assert parp1["evidence_type"] == "CRISPR/CRISPRi"

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_sl_pairs(str(tmp_path / "does_not_exist.tsv"))


class TestDGIdb:
    def test_approved_anticancer_drug_found(self, dgidb_tsv):
        idx = load_gene_to_drugs(str(dgidb_tsv))
        drugs = real_anticancer_drugs("PARP1", idx)
        names = {d["drug_name"] for d in drugs}
        assert names == {"OLAPARIB", "NIRAPARIB"}

    def test_unapproved_drug_excluded_by_default(self, dgidb_tsv):
        idx = load_gene_to_drugs(str(dgidb_tsv))
        drugs = real_anticancer_drugs("PARP1", idx)
        assert "EXPERIMENTAL-X" not in {d["drug_name"] for d in drugs}

    def test_corroboration_count_reflects_independent_sources(self, dgidb_tsv):
        idx = load_gene_to_drugs(str(dgidb_tsv))
        drugs = {d["drug_name"]: d for d in real_anticancer_drugs("PARP1", idx)}
        assert drugs["OLAPARIB"]["n_sources"] == 2  # ChEMBL + CIViC
        assert drugs["NIRAPARIB"]["n_sources"] == 1

    def test_biomarker_association_without_direct_action_excluded(self, dgidb_tsv):
        """Real regression case: a drug only ever reported as a bare
        biomarker/association claim (interaction_type NULL, no source ever
        calling it an inhibitor/antagonist/etc.) must not be returned --
        this is exactly what polluted the real CTLA4/PDCD1 lookups before
        the direct-action filter was added (see dgidb.py's module docstring)."""
        idx = load_gene_to_drugs(str(dgidb_tsv))
        drugs = real_anticancer_drugs("CTLA4", idx)
        names = {d["drug_name"] for d in drugs}
        assert "PREDNISONE" not in names  # NULL-type only, and not anti-neoplastic anyway
        assert "IPILIMUMAB" in names       # real inhibitor claim

    def test_require_direct_action_can_be_disabled(self, dgidb_tsv):
        idx = load_gene_to_drugs(str(dgidb_tsv))
        drugs = real_anticancer_drugs("PARP1", idx, require_direct_action=False)
        names = {d["drug_name"] for d in drugs}
        assert names == {"OLAPARIB", "NIRAPARIB"}  # both still direct-action here regardless

    def test_direct_action_types_is_lowercase(self):
        assert all(t == t.lower() for t in DIRECT_ACTION_TYPES)

    def test_real_immunotherapy_drugs(self, dgidb_tsv):
        out = real_immunotherapy_drugs(load_gene_to_drugs(str(dgidb_tsv)))
        names = {d["drug_name"] for d in out}
        assert names == {"IPILIMUMAB"}

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_gene_to_drugs(str(tmp_path / "nope.tsv"))


class TestAttachClinicalTrials:
    def _ct_file(self, tmp_path, results):
        path = tmp_path / "clinical_trials.json"
        path.write_text(json.dumps({"results": results}), encoding="utf-8")
        return str(path)

    def test_matching_row_gets_real_trials(self, tmp_path):
        trial = {"nct_id": "NCT00000001", "title": "Real Trial", "status": "RECRUITING"}
        ct = self._ct_file(tmp_path, {"hyperthermia_mild_41_42c::OLAPARIB": {"trials": [trial]}})
        rows = [{"modifier_id": "hyperthermia_mild_41_42c", "drug_name": "OLAPARIB"}]
        _attach_clinical_trials(rows, ct)
        assert rows[0]["has_clinical_trial"] is True
        assert rows[0]["clinical_trials"] == [trial]

    def test_unmatched_row_gets_empty_list_not_fabricated(self, tmp_path):
        ct = self._ct_file(tmp_path, {})
        rows = [{"modifier_id": "exercise", "drug_name": "NIVOLUMAB"}]
        _attach_clinical_trials(rows, ct)
        assert rows[0]["has_clinical_trial"] is False
        assert rows[0]["clinical_trials"] == []

    def test_missing_file_leaves_rows_untouched(self, tmp_path):
        rows = [{"modifier_id": "exercise", "drug_name": "NIVOLUMAB"}]
        _attach_clinical_trials(rows, str(tmp_path / "does_not_exist.json"))
        assert "clinical_trials" not in rows[0]
