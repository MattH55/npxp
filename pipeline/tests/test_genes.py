import numpy as np
import pytest

from npi_pharma.genes import GeneOverlapError, align, dedupe_vector, normalize_symbol


def test_normalize_maps_mouse_case_and_versions():
    assert normalize_symbol(" Mtor ") == "MTOR"
    assert normalize_symbol("ENSG00000198793.12") == "ENSG00000198793"


def test_dedupe_keeps_max_abs():
    genes, vals = dedupe_vector(["a", "A", "b"], np.array([0.5, -2.0, 1.0]))
    assert genes == ["A", "B"] and vals.tolist() == [-2.0, 1.0]


def test_alignment_raises_instead_of_going_empty():
    g1 = [f"G{i}" for i in range(500)]
    g2 = [f"H{i}" for i in range(500)]
    with pytest.raises(GeneOverlapError, match="0 genes"):
        align({"a": (g1, np.ones(500)), "b": (g2, np.ones(500))})


def test_alignment_threshold_boundary():
    a = [f"G{i}" for i in range(300)]
    b = [f"G{i}" for i in range(100, 400)]  # overlap 200
    shared, v = align({"a": (a, np.arange(300.0)), "b": (b, np.arange(300.0))})
    assert len(shared) == 200 and len(v["a"]) == 200
    with pytest.raises(GeneOverlapError):
        align({"a": (a, np.arange(300.0)), "b": (b[1:], np.arange(299.0))})


def test_alignment_preserves_values_by_gene():
    shared, v = align({"a": (["X", "Y"] + [f"G{i}" for i in range(200)], np.arange(202.0)),
                       "b": ([f"G{i}" for i in range(200)][::-1] + ["Y"], np.arange(201.0))}, min_overlap=10)
    i = shared.index("Y")
    assert v["a"][i] == 1.0 and v["b"][i] == 200.0
