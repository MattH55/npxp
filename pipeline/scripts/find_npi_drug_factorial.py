"""Are there GEO series with a factorial NPI x drug design? Measured: essentially none.

The models that predict drug synergy from expression (DeepSynergy, MatchMaker,
DRSPRING; docs/references.md) are all supervised on measured drug-PAIR synergy.
Applying anything of that shape to NPI x drug needs labels, and this script asks
whether public transcriptomics can supply them: a series profiling control, the NPI,
the drug, and the two together.

Result: of 78 series screened across hyperthermia+chemotherapy, fasting/caloric
restriction+chemotherapy, ketogenic/BHB+drug and hypoxia+drug queries, **2** had any
factorial design and **0** were NPI x drug. The one apparent hit, GSE290343, is a
false positive of the NPI term regex -- it matched "fast" inside "RNA-seq FASTQ were
mapped" in a processing field; its real arms are IL1B + 2-deoxyglucose, a cytokine
plus a drug (2DG is a caloric-restriction *mimetic*, which is NPI-adjacent
pharmacology, not an NPI).

How strong is that negative? Moderate, not conclusive. `detect_factorial` is strict --
it needs GEO's own annotations to separate control, single agents and combination
cleanly, which many series do not do -- and the queries are hand-written. So this
shows such designs are scarce in the obvious places, not that none exists. It is
enough to say that a supervised NPI x drug interaction model has no training set to
be fitted on today.

    python scripts/find_npi_drug_factorial.py
"""

import sys, json, time
sys.path.insert(0,'scripts'); sys.path.insert(0,'src')
from pathlib import Path
from find_combination_series import detect_factorial
from npi_pharma.ingest.fetch import fetch_geo_series
from npi_pharma.ingest.geo import read_series_matrix
import re
NPI=re.compile(r"hyperthermi|heat|42\s*°|43\s*°|thermal|fast|caloric|calorie|starv|"
               r"glucose|nutrient|serum.?free|deprivat|restrict|ketone|ketogenic|"
               r"hydroxybutyrate|hypoxia|anoxia", re.I)
data=json.load(open("out/npi_drug_factorial/candidates.json"))
targets=[]
for k in ["hyperthermia+chemo","fasting/CR+chemo","ketogenic/BHB+drug","hypoxia+drug"]:
    targets += [(k,r) for r in data[k]]
raw=Path("data/raw"); rows=[]
for k,r in targets:
    acc=r["accession"]
    try:
        fetch_geo_series(acc, raw, suppl=False, log=lambda *a,**kw: None)
        mats=sorted((raw/acc).glob("*series_matrix.txt.gz"))
        if not mats: continue
        _,s=read_series_matrix(mats[0])
    except Exception as e:
        rows.append({"q":k,"accession":acc,"outcome":f"error {type(e).__name__}"}); continue
    hit=detect_factorial(s)
    # does ANY annotation value name an NPI?
    npi_vals=set()
    for c in s.columns:
        if c in ("title","geo_accession"): continue
        for v in s[c].astype(str).unique():
            if NPI.search(v) and len(v)<90: npi_vals.add(f"{c}={v}")
    rows.append({"q":k,"accession":acc,"n":r["n"],
                 "factorial": bool(hit),
                 "field": hit["field"] if hit else None,
                 "singles": "; ".join(hit["singles"])[:90] if hit else None,
                 "combo": "; ".join(hit["combination"])[:70] if hit else None,
                 "npi_terms": "; ".join(sorted(npi_vals))[:130],
                 "title": r["title"][:90]})
    time.sleep(0.15)
import pandas as pd
t=pd.DataFrame(rows)
t.to_csv("out/npi_drug_factorial/screened.tsv",sep="\t",index=False)
print(f"screened {len(t)}")
f=t[t.get("factorial")==True]
print(f"factorial designs: {len(f)}")
g=f[f.npi_terms.astype(str).str.len()>0]
print(f"factorial AND naming an NPI: {len(g)}\n")
for r in g.itertuples():
    print(f"== {r.accession} [{r.q}] n={r.n}\n   field: {r.field}\n   singles: {r.singles}\n   combo: {r.combo}\n   npi: {r.npi_terms}\n   {r.title}\n")
