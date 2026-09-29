"""Run the drug-consensus build one drug per process, then merge the parts.

`build_drug_consensus.py` holds every drug's per-series signatures in memory so it
can compute cross-series agreement at the end. With the RNA-seq path enabled that
is 20 drugs x up to 10 series x ~30,000 genes, and the run was killed by the OOM
killer partway through carboplatin -- silently, with no traceback, which is what an
out-of-memory kill looks like from inside the log.

Nothing about the method needs the drugs to share a process: a consensus is taken
within a drug, and the cross-series agreement that feeds the reliability bands is
computed within a drug too. So each drug runs on its own and the parts are merged
here. Memory is then bounded by the largest single drug rather than the sum, and a
drug that dies takes only itself down -- the rest of the run survives, and rerunning
skips the parts already on disk.

    python scripts/build_drug_consensus_all.py
    python scripts/build_drug_consensus_all.py --drugs cisplatin oxaliplatin --force
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from build_drug_consensus import MISSING_FROM_LINCS  # noqa: E402

from npi_pharma.store import load_signatures, save_signatures  # noqa: E402


def slug(drug: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", drug.lower()).strip("_")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--drugs", nargs="*", default=MISSING_FROM_LINCS)
    ap.add_argument("--max-per-drug", type=int, default=80)
    ap.add_argument("--max-series-used", type=int, default=10)
    ap.add_argument("--raw-dir", default="data/raw")
    ap.add_argument("--out", default="data/processed/signatures/drugs_consensus_v2.parquet")
    ap.add_argument("--report", default="out/drug_consensus_v2")
    ap.add_argument("--timeout", type=int, default=3600, help="seconds per drug")
    ap.add_argument("--min-free-gb", type=float, default=12.0,
                    help="stop before a drug that would run with less headroom than this; "
                         "the supplementary downloads are several GB per drug, and running "
                         "out mid-write loses the part rather than skipping it")
    ap.add_argument("--force", action="store_true", help="rebuild parts already on disk")
    a = ap.parse_args(argv)
    rep = Path(a.report)
    parts = rep / "parts"
    parts.mkdir(parents=True, exist_ok=True)
    here = Path(__file__).parent

    status = []
    for drug in a.drugs:
        d = parts / slug(drug)
        done = d / "consensus.parquet"
        if done.exists() and not a.force:
            status.append({"drug": drug, "outcome": "kept existing part"})
            print(f"== {drug}: part already built, kept", file=sys.stderr)
            continue
        free_gb = shutil.disk_usage(a.raw_dir).free / 2**30
        if free_gb < a.min_free_gb:
            status.append({"drug": drug, "outcome": f"skipped: only {free_gb:.1f} GB free"})
            print(f"== {drug}: only {free_gb:.1f} GB free, stopping here", file=sys.stderr)
            break
        d.mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, str(here / "build_drug_consensus.py"), "--allow-rnaseq",
               "--drugs", drug, "--max-per-drug", str(a.max_per_drug),
               "--max-series-used", str(a.max_series_used), "--raw-dir", a.raw_dir,
               "--out", str(done), "--report", str(d)]
        t0 = time.time()
        print(f"== {drug}", file=sys.stderr)
        try:
            p = subprocess.run(cmd, timeout=a.timeout, capture_output=True, text=True)
            (d / "run.log").write_text(p.stdout + "\n" + p.stderr)
            # A negative return code is a signal: -9 is the OOM killer, which is the
            # failure this script exists to contain. Record it rather than hide it.
            outcome = ("built" if p.returncode == 0 else
                       f"killed by signal {-p.returncode}" if p.returncode < 0 else
                       f"exit {p.returncode}")
        except subprocess.TimeoutExpired:
            outcome = f"timed out after {a.timeout}s"
        status.append({"drug": drug, "outcome": outcome,
                       "seconds": round(time.time() - t0, 1),
                       "part_built": done.exists()})
        print(f"   {outcome} in {time.time() - t0:.0f}s", file=sys.stderr)

    # merge: consensuses, per-series signatures, series logs and QC
    cons, per_series, logs, qc = [], [], [], {}
    for drug in a.drugs:
        d = parts / slug(drug)
        if (d / "consensus.parquet").exists():
            cons += list(load_signatures(d / "consensus.parquet"))
        ps = d / "drugs_geo_per_series.parquet"
        if ps.exists():
            per_series += list(load_signatures(ps))
        sl = d / "series_log.tsv"
        if sl.exists() and sl.stat().st_size:
            logs.append(pd.read_csv(sl, sep="\t"))
        q = d / "qc.json"
        if q.exists():
            qc |= json.loads(q.read_text())

    if cons:
        save_signatures(cons, a.out)
    if per_series:
        save_signatures(per_series, str(Path(a.out).with_name("drugs_geo_per_series_v2.parquet")))
    if logs:
        pd.concat(logs, ignore_index=True).to_csv(rep / "series_log.tsv", sep="\t", index=False)
    (rep / "qc.json").write_text(json.dumps(qc, indent=2, default=float))
    pd.DataFrame(status).to_csv(rep / "run_status.tsv", sep="\t", index=False)

    print(f"\n{len(cons)} consensus signatures from {len(per_series)} series -> {a.out}")
    if cons:
        rows = [{"drug": c.sig_id, "n_series": (c.meta or {}).get("n_series"),
                 "n_genes": len(c.genes),
                 "cross_series": (qc.get(c.sig_id) or {}).get("median_cross_series_cosine"),
                 "vs_lincs": (qc.get(c.sig_id) or {}).get("vs_lincs_consensus_cosine")}
                for c in cons]
        t = pd.DataFrame(rows).sort_values("n_series", ascending=False)
        print(t.to_string(index=False, na_rep="-"))
    bad = [s for s in status if s.get("outcome") not in ("built", "kept existing part")]
    if bad:
        print("\ndrugs that did not complete:")
        print(pd.DataFrame(bad).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
