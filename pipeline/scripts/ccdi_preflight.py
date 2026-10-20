"""Probe the CCDI data hosts named in docs/ccdi_integration_brief.md section 1.

    python scripts/ccdi_preflight.py

Prints reachable / blocked per host. Blocked hosts are reported, never worked
around with mirrors of unknown provenance.
"""

from __future__ import annotations

import urllib.error
import urllib.request

HOSTS = {
    "GDC (TARGET, HCMI)": "https://api.gdc.cancer.gov/status",
    "refine.bio": "https://api.refine.bio/v1/",
    "Treehouse": "https://treehousegenomics.soe.ucsc.edu/public-data/",
    "OpenPedCan bucket": "https://s3.amazonaws.com/d3b-openaccess-us-east-1-prd-pbta/open-targets/v15/md5sum.txt",
    "Molecular Targets Platform": "https://moleculartargets.ccdi.cancer.gov/",
    "ScPCA": "https://api.scpca.alexslemonade.org/v1/projects/",
    "PedcBioPortal": "https://pedcbioportal.kidsfirstdrc.org/api/health",
    "FDA (PMTL)": "https://www.fda.gov/",
    "SEER": "https://seer.cancer.gov/",
    "NCBI GEO": "https://ftp.ncbi.nlm.nih.gov/geo/",
}


def probe(url: str, timeout: int = 15) -> str:
    try:
        with urllib.request.urlopen(urllib.request.Request(url, method="GET"), timeout=timeout) as r:
            return f"reachable ({r.status})"
    except urllib.error.HTTPError as e:  # the host answered
        return f"reachable (HTTP {e.code})"
    except Exception as e:  # proxy denial, DNS, timeout
        return f"BLOCKED ({getattr(e, 'reason', e)})"


if __name__ == "__main__":
    for name, url in HOSTS.items():
        print(f"{name:28s} {probe(url)}")
