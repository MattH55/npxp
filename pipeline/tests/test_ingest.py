import gzip

import h5py
import numpy as np
import pandas as pd
import pytest

from npi_pharma.ingest.geo import collapse_probes, read_series_matrix, split_pre_post
from npi_pharma.ingest.lincs import ingest_lincs, read_gctx_columns


def _write_gctx(path, genes, sigs, mat):
    """Minimal GCTX (HDF5) in the CMap layout: matrix stored as [col, row]."""
    with h5py.File(path, "w") as f:
        f["0/DATA/0/matrix"] = mat.T.astype(np.float32)
        f["0/META/ROW/id"] = np.array(genes, dtype="S")
        f["0/META/COL/id"] = np.array(sigs, dtype="S")


@pytest.fixture
def lincs_files(tmp_path):
    rng = np.random.default_rng(0)
    entrez = [str(1000 + i) for i in range(400)]
    sigs = [f"SIG{i}" for i in range(12)]
    mat = rng.normal(size=(400, 12))
    mat[0, :4] = 5.0  # metformin signatures share a strong gene
    _write_gctx(tmp_path / "l5.gctx", entrez, sigs, mat)
    si = pd.DataFrame({
        "sig_id": sigs,
        "pert_id": ["BRD-MET"] * 4 + ["BRD-SIR"] * 4 + [f"BRD-{i}" for i in range(4)],
        "pert_iname": ["metformin"] * 4 + ["sirolimus"] * 4 + [f"cmpd{i}" for i in range(4)],
        "pert_type": ["trt_cp"] * 12,
        "cell_id": ["MCF7", "A549", "PC3", "HEPG2"] * 3,
        "pert_idose": ["10 uM"] * 12,
        "pert_itime": ["24 h"] * 12,
    })
    si.to_csv(tmp_path / "sig_info.txt", sep="\t", index=False)
    gi = pd.DataFrame({"pr_gene_id": entrez, "pr_gene_symbol": [f"GENE{i}" for i in range(400)],
                       "pr_is_lm": ["1"] * 100 + ["0"] * 300, "pr_is_bing": ["1"] * 300 + ["0"] * 100})
    gi.to_csv(tmp_path / "gene_info.txt", sep="\t", index=False)
    return tmp_path, mat


def test_read_gctx_columns_orders_and_values(lincs_files):
    d, mat = lincs_files
    df = read_gctx_columns(d / "l5.gctx", ["SIG5", "SIG1"])
    assert list(df.columns) == ["SIG1", "SIG5"]
    np.testing.assert_allclose(df["SIG5"].to_numpy(), mat[:, 5], rtol=1e-6)


def test_ingest_lincs_consensus_alias_and_missing(lincs_files):
    d, mat = lincs_files
    sigs, missing = ingest_lincs(d / "l5.gctx", d / "sig_info.txt", d / "gene_info.txt",
                                 ["metformin", "rapamycin", "notadrug"], gene_space="bing")
    by = {s.sig_id: s for s in sigs}
    assert set(by) == {"metformin", "rapamycin"} and missing == ["notadrug"]
    met = by["metformin"]
    assert len(met.genes) == 300 and met.meta["n_signatures"] == 4
    assert met.z[met.genes.index("GENE0")] == pytest.approx(5.0)
    np.testing.assert_allclose(met.z[met.genes.index("GENE7")], np.median(mat[7, :4]), rtol=1e-5)
    assert by["rapamycin"].meta["pert_iname"] == "sirolimus"


def test_ingest_lincs_cell_line_filter_and_random(lincs_files):
    d, _ = lincs_files
    sigs, missing = ingest_lincs(d / "l5.gctx", d / "sig_info.txt", d / "gene_info.txt", ["metformin"],
                                 cell_lines=["MCF7"], gene_space="landmark", n_random=2, seed=1)
    met = next(s for s in sigs if s.sig_id == "metformin")
    assert met.tissue == "MCF7" and len(met.genes) == 100
    # random picks with no MCF7 signature are reported as missing, not faked
    assert len(sigs) + len(missing) == 3 and "metformin" not in missing
    again, _ = ingest_lincs(d / "l5.gctx", d / "sig_info.txt", d / "gene_info.txt", ["metformin"],
                            cell_lines=["MCF7"], gene_space="landmark", n_random=2, seed=1)
    assert [s.sig_id for s in sigs] == [s.sig_id for s in again]


SERIES = """!Series_title\t"toy"
!Sample_title\t"p1 pre"\t"p1 post"\t"p2 pre"\t"p2 post"
!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"\t"GSM4"
!Sample_characteristics_ch1\t"subject: p1"\t"subject: p1"\t"subject: p2"\t"subject: p2"
!Sample_characteristics_ch1\t"time: baseline"\t"time: week 8"\t"time: baseline"\t"time: week 8"
!series_matrix_table_begin
"ID_REF"\t"GSM1"\t"GSM2"\t"GSM3"\t"GSM4"
"pr1"\t5\t6\t5.5\t6.5
"pr2"\t7\t7\t7\t7
"pr3"\t1\t2\t1\t2
!series_matrix_table_end
"""


def test_series_matrix_pairing(tmp_path):
    p = tmp_path / "toy_series_matrix.txt.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(SERIES)
    expr, samples = read_series_matrix(p)
    assert expr.shape == (3, 4) and list(samples["time"]) == ["baseline", "week 8"] * 2
    genes = collapse_probes(expr, pd.Series({"pr1": "MTOR", "pr2": "RPTOR", "pr3": "MTOR"}))
    assert sorted(genes.index) == ["MTOR", "RPTOR"]
    assert genes.loc["MTOR", "GSM1"] == 5  # highest-mean probe kept
    pre, post, paired = split_pre_post(genes, samples, {"subject_field": "subject", "time_field": "time",
                                                        "pre": "baseline", "post": "week 8"})
    assert paired and list(pre.columns) == ["p1", "p2"] and post.loc["MTOR", "p2"] == 6.5
