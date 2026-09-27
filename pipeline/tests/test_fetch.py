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
