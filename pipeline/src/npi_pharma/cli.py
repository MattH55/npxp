"""npi-pharma command-line interface."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import yaml

from .config import load_config
from .eval.report import envelope, explain, write_json, write_tsv
from .gene_sets import load_gmt
from .interact.rank import rank_pairs, select_npis
from .store import load_patient, load_signatures, save_patient, save_signatures, signature_index


def _drug_set(path: str | None, name: str | None) -> list[str] | None:
    if not name:
        return None
    sets = yaml.safe_load(Path(path).read_text())["drug_sets"]
    if name not in sets:
        raise SystemExit(f"drug set {name!r} not in {path}; have {sorted(sets)}")
    return [d.lower() for d in sets[name]]


def cmd_ingest_lincs(a) -> int:
    from .ingest.lincs import ingest_lincs

    if a.level != 5:
        raise SystemExit("only Level 5 (MODZ consensus) signatures are supported")
    perts = [p.strip() for p in a.perts.split(",")] if a.perts else _drug_set(a.drug_sets, a.drug_set)
    if not perts:
        raise SystemExit("pass --perts or --drug-set")
    targets = None
    if a.targets:
        t = pd.read_csv(a.targets, sep="\t", dtype=str)
        targets = {r.pert_iname.lower(): str(r.target).split("|") for r in t.itertuples()}
    sigs, missing = ingest_lincs(
        a.gctx, a.sig_info, a.gene_info, perts, sig_metrics=a.sig_metrics, pert_info=a.pert_info,
        targets=targets, cell_lines=a.cell_lines.split(",") if a.cell_lines else None,
        time_h=[int(t) for t in a.time_h.split(",")] if a.time_h else None,
        min_tas=a.min_tas, gene_space=a.gene_space, n_random=a.random, seed=a.seed,
    )
    save_signatures(sigs, a.out)
    print(f"wrote {len(sigs)} drug signatures -> {a.out}")
    if missing:
        print(f"not found in sig_info ({len(missing)}): {', '.join(missing)}", file=sys.stderr)
    return 0


def cmd_ingest_npi(a) -> int:
    from .ingest.catalog import build_from_entry, load_catalog

    cfg = load_config(a.config)
    entries = load_catalog(a.catalog, min_n=cfg["min_npi_n"])
    ids = set(a.ids.split(",")) if a.ids else None
    built, skipped = [], []
    for e in entries:
        if ids and e.npi_id not in ids:
            continue
        try:
            built.append(build_from_entry(e, a.raw_dir))
        except (ValueError, FileNotFoundError, KeyError) as err:
            skipped.append(str(err) if str(err).startswith(e.npi_id) else f"{e.npi_id}: {err}")
    if built:
        save_signatures(built, a.out)
    print(f"built {len(built)} NPI signatures" + (f" -> {a.out}" if built else ""))
    for s in skipped:
        print(f"  skipped {s}", file=sys.stderr)
    for e in entries:
        if e.issues and (not ids or e.npi_id in ids):
            print(f"  catalog {e.npi_id}: {'; '.join(e.issues)}", file=sys.stderr)
    return 0 if built or not skipped else 1


def cmd_encode_patient(a) -> int:
    from .patient.encode import encode_patient, read_expression

    obs_filter = dict(kv.split("=", 1) for kv in a.obs_filter) if a.obs_filter else None
    expr, flags = read_expression(a.input, obs_filter=obs_filter, sample_key=a.sample_key)
    ref = read_expression(a.reference)[0] if a.reference else None
    p = encode_patient(expr, a.sample, a.tissue, reference=ref, allow_no_reference=a.allow_no_reference, flags=flags)
    save_patient(p, a.out)
    print(f"encoded {p.sample_id} ({p.tissue}, {len(p.genes)} genes; s_P vs {p.reference}) -> {a.out}")
    for f in p.flags:
        print(f"  flag: {f}", file=sys.stderr)
    return 0


def _load_common(a):
    cfg = load_config(a.config)
    return cfg, load_gmt(a.gene_sets), load_patient(a.patient), load_signatures(a.npis), load_signatures(a.drugs)


def _inputs(a, **extra) -> dict:
    return {"patient": a.patient, "npis": a.npis, "drugs": a.drugs, "gene_sets": a.gene_sets or "bundled curated_core.gmt",
            "config": a.config or "bundled scoring.yaml", **extra}


def cmd_score(a) -> int:
    cfg, gs, patient, npis, drugs = _load_common(a)
    nidx, didx = signature_index(npis), signature_index(drugs)
    if a.npi not in nidx:
        raise SystemExit(f"NPI {a.npi!r} not in {a.npis}")
    if a.drug not in didx:
        raise SystemExit(f"drug {a.drug!r} not in {a.drugs}")
    # Rank the query against every drug for this NPI so composite_rank means something.
    table, reports = rank_pairs(patient, [nidx[a.npi]], drugs, gs, cfg)
    rep = reports.get((a.npi, a.drug))
    if rep is None:
        row = table[(table.npi_id == a.npi) & (table.drug_id == a.drug)]
        raise SystemExit(f"could not score pair: {row['flags'].iat[0] if len(row) else 'unknown'}")
    top = table.dropna(subset=["composite"]).iloc[0]
    body = {
        "patient": {"sample_id": patient.sample_id, "tissue": patient.tissue, "reference": patient.reference,
                    "flags": patient.flags},
        "score": rep,
        "explanation": explain(rep),
        "top_pair_for_this_npi": explain(reports[(top.npi_id, top.drug_id)]),
    }
    write_json(envelope(body, _inputs(a, npi=a.npi, drug=a.drug)), a.out)
    if a.tsv:
        write_tsv(table, a.tsv)
    print(f"{a.npi} x {a.drug}: composite {rep['composite']:+.4f} (rank {rep['composite_rank']}/{rep['n_pairs_ranked']}), "
          f"confidence x{rep['confidence']:.2f} -> {a.out}")
    for w in rep["warnings"]:
        print(f"  warning: {w}", file=sys.stderr)
    return 0


def cmd_rank(a) -> int:
    cfg, gs, patient, npis, drugs = _load_common(a)
    ids = a.npi.split(",") if a.npi else None
    sel, skipped = select_npis(npis, patient, cfg, npi_class=a.npi_class, npi_ids=ids, allow_mismatch=a.allow_mismatch)
    for s in skipped:
        print(f"  skipped {s}", file=sys.stderr)
    if not sel:
        raise SystemExit("no NPI signatures selected (check --npi-class / tissue compatibility / --allow-mismatch)")
    wanted = _drug_set(a.drug_sets, a.drug_set)
    if wanted:
        drugs = [d for d in drugs if d.sig_id.lower() in wanted or d.meta.get("pert_iname", "").lower() in wanted]
        if not drugs:
            raise SystemExit(f"no drugs in {a.drugs} match drug set {a.drug_set!r}")
    table, reports = rank_pairs(patient, sel, drugs, gs, cfg)
    write_tsv(table.head(a.top) if a.top else table, a.out)
    if a.json:
        scored = table.dropna(subset=["composite"])
        top = [explain(reports[(r.npi_id, r.drug_id)]) for r in scored.head(min(a.top or 5, 5)).itertuples()]
        write_json(envelope({"n_pairs": len(table), "top_pairs": top},
                            _inputs(a, npi_class=a.npi_class, npi=a.npi, drug_set=a.drug_set)), a.json)
    with pd.option_context("display.width", 200, "display.max_columns", 12):
        cols = ["composite_rank", "npi_id", "drug_id", "composite", "complementarity", "pathway_joint", "confidence"]
        print(table[cols].head(a.top or 25).to_string(index=False))
    return 0


def cmd_make_fixture(a) -> int:
    from .fixtures import write_fixture

    for k, p in write_fixture(a.out).items():
        print(f"{k}: {p}")
    print("NOTE: synthetic fixture data (FIXTURE_* ids), for tests and demos only.")
    return 0


def cmd_demo(a) -> int:
    """Offline end-to-end run of the milestone path on synthetic fixtures."""
    out = Path(a.out)
    fx = dict(npis=out / "npis.parquet", drugs=out / "drugs.parquet", patient=out / "patient_adipose.tsv",
              reference=out / "reference_adipose.tsv")
    cmd_make_fixture(argparse.Namespace(out=out))
    cmd_encode_patient(argparse.Namespace(input=fx["patient"], reference=fx["reference"], sample=None, tissue="adipose",
                                          allow_no_reference=False, out=out / "patient.npz", obs_filter=None, sample_key=None))
    common = dict(patient=out / "patient.npz", npis=fx["npis"], drugs=fx["drugs"], config=None, gene_sets=None)
    cmd_score(argparse.Namespace(**common, npi="FIXTURE_LCD_adipose", drug="FIXTURE_metformin",
                                 out=out / "report.json", tsv=out / "ranked.tsv"))
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="npi-pharma", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("ingest-lincs", help="LINCS L1000 Level 5 GCTX -> drug signature parquet")
    s.add_argument("--level", type=int, default=5)
    s.add_argument("--gctx", required=True)
    s.add_argument("--sig-info", required=True)
    s.add_argument("--gene-info", required=True)
    s.add_argument("--sig-metrics")
    s.add_argument("--pert-info")
    s.add_argument("--targets", help="TSV with pert_iname, target (|-separated), e.g. from the Repurposing Hub")
    s.add_argument("--perts", help="comma-separated pert_iname list")
    s.add_argument("--drug-set")
    s.add_argument("--drug-sets", default="configs/drug_sets.yaml")
    s.add_argument("--cell-lines")
    s.add_argument("--time-h")
    s.add_argument("--min-tas", type=float)
    s.add_argument("--gene-space", choices=["landmark", "bing", "all"], default="bing")
    s.add_argument("--random", type=int, default=0, help="add N random trt_cp compounds (background)")
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_ingest_lincs)

    s = sub.add_parser("ingest-npi", help="build NPI signatures from catalog entries with local inputs")
    s.add_argument("--catalog", default="configs/npi_catalog.yaml")
    s.add_argument("--raw-dir", default="data/raw")
    s.add_argument("--ids")
    s.add_argument("--config")
    s.add_argument("--out", default="data/processed/signatures/npis.parquet")
    s.set_defaults(func=cmd_ingest_npi)

    s = sub.add_parser("encode-patient", help="expression -> PatientState (.npz)")
    s.add_argument("--input", required=True, help="genes x samples TSV/CSV, or .h5ad")
    s.add_argument("--sample")
    s.add_argument("--tissue", required=True)
    s.add_argument("--reference", help="healthy reference, genes x samples (same tissue/platform)")
    s.add_argument("--allow-no-reference", action="store_true")
    s.add_argument("--sample-key", help=".h5ad obs column to pseudobulk by")
    s.add_argument("--obs-filter", nargs="*", help=".h5ad obs filters key=value")
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_encode_patient)

    for name, fn, hlp in (("score", cmd_score, "score one NPI x drug pair"), ("rank", cmd_rank, "rank NPI x drug pairs")):
        s = sub.add_parser(name, help=hlp)
        s.add_argument("--patient", required=True)
        s.add_argument("--npis", default="data/processed/signatures/npis.parquet")
        s.add_argument("--drugs", default="data/processed/signatures/drugs.parquet")
        s.add_argument("--gene-sets", help="GMT file (default: bundled curated core sets)")
        s.add_argument("--config", help="scoring YAML overriding bundled defaults")
        s.add_argument("--out", required=True)
        s.set_defaults(func=fn)
    sp = sub.choices["score"]
    sp.add_argument("--npi", required=True)
    sp.add_argument("--drug", required=True)
    sp.add_argument("--tsv", help="also write the ranked table for this NPI")
    rp = sub.choices["rank"]
    rp.add_argument("--npi", help="comma-separated NPI ids (bypasses tissue filter)")
    rp.add_argument("--npi-class", help="NPI modality, e.g. exercise, diet, caloric_restriction")
    rp.add_argument("--drug-set")
    rp.add_argument("--drug-sets", default="configs/drug_sets.yaml")
    rp.add_argument("--top", type=int, default=25)
    rp.add_argument("--allow-mismatch", action="store_true", help="keep tissue-mismatched NPIs (down-weighted)")
    rp.add_argument("--json", help="also write JSON explanations for the top pairs")

    s = sub.add_parser("make-fixture", help="write the synthetic 50-drug fixture")
    s.add_argument("--out", default="data/fixture")
    s.set_defaults(func=cmd_make_fixture)
    s = sub.add_parser("demo", help="offline milestone path on synthetic fixtures")
    s.add_argument("--out", default="out/demo")
    s.set_defaults(func=cmd_demo)
    return ap


def main(argv: list[str] | None = None) -> int:
    a = build_parser().parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    raise SystemExit(main())
