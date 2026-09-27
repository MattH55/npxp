import pytest

from npi_pharma.config import load_config
from npi_pharma.fixtures import build_all
from npi_pharma.gene_sets import load_gmt
from npi_pharma.patient.encode import encode_patient


@pytest.fixture(scope="session")
def fx():
    return build_all()


@pytest.fixture(scope="session")
def cfg():
    return load_config()


@pytest.fixture(scope="session")
def gene_sets():
    return load_gmt()


@pytest.fixture(scope="session")
def patient(fx):
    return encode_patient(fx["patient_expr"], None, "adipose", reference=fx["reference_expr"])


@pytest.fixture(scope="session")
def by_id(fx):
    return {s.sig_id: s for s in fx["npis"] + fx["drugs"]}
