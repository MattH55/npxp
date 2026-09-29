import gzip

import pytest

from npi_pharma.ingest import fetch


def test_geo_directory_stems():
    assert fetch.series_dir("GSE95640", "matrix").endswith("/series/GSE95nnn/GSE95640/matrix/")
    assert fetch.series_dir("GSE770", "suppl").endswith("/series/GSEnnn/GSE770/suppl/")
    assert fetch.platform_annot_url("GPL6244").endswith("/platforms/GPL6nnn/GPL6244/annot/GPL6244.annot.gz")
    assert "/platforms/GPLnnn/GPL570/" in fetch.platform_annot_url("GPL570")


def test_annot_to_probe_map(tmp_path):
    annot = tmp_path / "GPL1.annot.gz"
    with gzip.open(annot, "wt") as fh:
        fh.write("^Annotation\n!platform_table_begin\nID\tIDENTIFIER\tGene symbol\n"
                 "p1\tx\tMTOR\np2\tx\t\np3\tx\tRPTOR\n!platform_table_end\n")
    out = fetch.annot_to_probe_map(annot, tmp_path / "pm.tsv")
    assert out.read_text().splitlines() == ["p1\tMTOR", "p3\tRPTOR"]


def test_fetch_geo_series_offline(tmp_path, monkeypatch):
    matrix = gzip.compress(b'!Series_platform_id\t"GPL1"\n!series_matrix_table_begin\n')
    annot = gzip.compress(b"!platform_table_begin\nID\tGene symbol\np1\tMTOR\n!platform_table_end\n")
    pages = {
        fetch.series_dir("GSE1234", "matrix"): b'<a href="GSE1234_series_matrix.txt.gz">x</a><a href="/geo/">up</a>',
        fetch.series_dir("GSE1234", "matrix") + "GSE1234_series_matrix.txt.gz": matrix,
        fetch.platform_annot_url("GPL1"): annot,
    }

    class Resp:
        def __init__(self, b):
            import io
            self.b = io.BytesIO(b)

        def read(self, *a):
            return self.b.read(*a)

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    def fake_urlopen(url, timeout=None):
        if url not in pages:
            raise OSError(f"404 {url}")
        return Resp(pages[url])

    monkeypatch.setattr(fetch.urllib.request, "urlopen", fake_urlopen)
    got = fetch.fetch_geo_series("GSE1234", tmp_path, log=lambda *a, **k: None)
    assert got["series_matrix"][0].name == "GSE1234_series_matrix.txt.gz"
    assert (tmp_path / "GSE1234" / "probe_map.tsv").read_text() == "p1\tMTOR\n"
    with pytest.raises(OSError):  # suppl listing missing
        fetch.fetch_lincs("GSE1234", tmp_path, log=lambda *a, **k: None)


def test_gene_info_to_id_map_prefers_protein_coding(tmp_path):
    gi = tmp_path / "gi.gz"
    with gzip.open(gi, "wt") as fh:
        fh.write("#tax_id\tGeneID\tSymbol\tdbXrefs\ttype_of_gene\n"
                 "9606\t1\tMTOR-AS\tEnsembl:ENSG1\tncRNA\n"
                 "9606\t2\tMTOR\tMIM:1|HGNC:HGNC:2|Ensembl:ENSG1\tprotein-coding\n"
                 "9606\t3\tRPTOR\tEnsembl:ENSG2|Ensembl:ENSG3\tprotein-coding\n"
                 "9606\t4\tNOENS\t-\tprotein-coding\n")
    out = fetch.gene_info_to_id_map(gi, tmp_path / "m.tsv")
    assert out.read_text().splitlines() == ["ENSG1\tMTOR", "ENSG2\tRPTOR", "ENSG3\tRPTOR"]


def test_is_suppl_table_keeps_tables_and_rejects_unparseable_bulk():
    """The filter that stopped one series from filling the disk.

    GSE236253 ships a 2.7 GB Hi-C contact map beside its 4.8 MB count table, and
    GSE276609 ships 8.5 GB of 10x Loupe projects. Neither can be read by
    supplementary_table(), so neither should ever be fetched.
    """
    from npi_pharma.ingest.fetch import is_suppl_table

    for name in ["GSE1_counts.txt.gz", "GSE1_Read_counts.csv.gz", "GSE1_gene_fpkm.tsv",
                 "GSE1_table.xlsx", "GSE1_expression.txt", "GSE1_tpm.csv",
                 "GSE1_matrix.mtx.gz", "GSE1_counts.tab.bz2"]:
        assert is_suppl_table(name), name

    for name in ["GSE236253_U87.allValidPairs.hic", "GSE276609_Cancer.cloupe.gz",
                 "GSE276609_reanalysis.tar.gz", "GSE1_RAW.tar", "GSE1.bam", "GSE1.bw",
                 "GSE1_signal.bigWig", "GSE1.h5", "GSE1_img.png", "GSE1.pdf",
                 "GSE1_fastq.gz", "GSE1.CEL.gz", "GSE1.idat",
                 "GSE1_barcodes.tsv.tar.gz"]:   # a table inside a tar is not readable either
        assert not is_suppl_table(name), name


def test_download_abandons_a_file_over_its_cap(tmp_path, monkeypatch):
    """A cap is enforced on the declared length AND on the stream itself."""
    import io
    import urllib.request

    import pytest

    from npi_pharma.ingest import fetch

    class Resp(io.BytesIO):
        def __init__(self, data, declared):
            super().__init__(data)
            self.headers = {"Content-Length": str(declared)} if declared is not None else {}

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    payload = b"x" * (3 << 20)          # 3 MB

    # declared too large -> refused before a byte is written
    monkeypatch.setattr(urllib.request, "urlopen",
                        lambda *a, **k: Resp(payload, 3 << 20))
    dest = tmp_path / "big.txt.gz"
    with pytest.raises(fetch.TooLarge):
        fetch.download("http://x/big.txt.gz", dest, max_bytes=1 << 20)
    assert not dest.exists() and not dest.with_suffix(".gz.part").exists()

    # no declared length -> the stream itself is capped, and the partial file removed
    monkeypatch.setattr(urllib.request, "urlopen", lambda *a, **k: Resp(payload, None))
    with pytest.raises(fetch.TooLarge):
        fetch.download("http://x/big.txt.gz", dest, max_bytes=1 << 20)
    assert not dest.exists() and not dest.with_suffix(".gz.part").exists()

    # within the cap -> written
    monkeypatch.setattr(urllib.request, "urlopen", lambda *a, **k: Resp(payload, None))
    fetch.download("http://x/ok.txt.gz", tmp_path / "ok.txt.gz", max_bytes=8 << 20)
    assert (tmp_path / "ok.txt.gz").stat().st_size == len(payload)
