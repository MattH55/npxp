"""LINCS L1000 Level 5 (MODZ consensus z-score) ingestion from GCTX.

Works with the GEO-hosted releases (no clue.io login needed):
  Phase I  GSE92742: GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx
  Phase II GSE70138: GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_<date>.gctx
plus the matching ``*_sig_info*.txt``, ``*_gene_info*.txt`` and optionally
``*_sig_metrics*.txt`` (TAS) and ``*_pert_info*.txt`` (SMILES).

Only the columns for the requested perturbagens are read from the HDF5 file,
so a curated drug set costs megabytes of RAM, not the full 5-20 GB matrix.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd

from ..genes import dedupe_vector
from ..model import DRUG, PROV_LINCS, Signature

GCTX_MATRIX = "0/DATA/0/matrix"
GCTX_ROW_IDS = "0/META/ROW/id"
GCTX_COL_IDS = "0/META/COL/id"

# CMap names some compounds differently from common usage.
DEFAULT_ALIASES = {"rapamycin": ["sirolimus"], "sirolimus": ["rapamycin"]}


def _decode(arr) -> list[str]:
    return [a.decode() if isinstance(a, bytes) else str(a) for a in arr]


def read_gctx_columns(gctx: str | Path, col_ids: Iterable[str]) -> pd.DataFrame:
    """Return a genes x signatures DataFrame for the requested signature ids."""
    import h5py

    wanted = list(dict.fromkeys(col_ids))
    with h5py.File(gctx, "r") as f:
        cols = _decode(f[GCTX_COL_IDS][:])
        rows = _decode(f[GCTX_ROW_IDS][:])
        pos = {c: i for i, c in enumerate(cols)}
        missing = [c for c in wanted if c not in pos]
        if missing:
            raise KeyError(f"{len(missing)} signature ids not in GCTX, e.g. {missing[:3]}")
        idx = sorted(pos[c] for c in wanted)  # h5py fancy indexing needs increasing order
        # GCTX stores the matrix transposed: dataset[col, row].
        block = f[GCTX_MATRIX][idx, :]
    return pd.DataFrame(block.T, index=rows, columns=[cols[i] for i in idx])


def load_gene_map(gene_info: str | Path, space: str = "bing") -> dict[str, str]:
    """Entrez id -> symbol. ``space``: landmark (978), bing (landmark + best
    inferred, ~10k; the brief's "LINCS inferred genes"), or all (12,328)."""
    gi = pd.read_csv(gene_info, sep="\t", dtype=str)
    if space == "landmark":
        gi = gi[gi["pr_is_lm"].astype(str) == "1"]
    elif space == "bing":
        gi = gi[gi["pr_is_bing"].astype(str) == "1"]
    elif space != "all":
        raise ValueError(f"unknown gene space {space!r}")
    return dict(zip(gi["pr_gene_id"], gi["pr_gene_symbol"]))


def select_signatures(
    sig_info: pd.DataFrame,
    perts: Iterable[str],
    cell_lines: Iterable[str] | None = None,
    time_h: Iterable[int] | None = None,
    min_tas: float | None = None,
    aliases: dict[str, list[str]] | None = None,
) -> dict[str, pd.DataFrame]:
    """Map each requested compound name -> the sig_info rows to collapse."""
    aliases = DEFAULT_ALIASES if aliases is None else aliases
    si = sig_info[sig_info["pert_type"] == "trt_cp"] if "pert_type" in sig_info else sig_info
    names = si["pert_iname"].str.lower()
    out: dict[str, pd.DataFrame] = {}
    for p in perts:
        keys = [p.lower()] + [a.lower() for a in aliases.get(p.lower(), [])]
        rows = si[names.isin(keys)]
        if cell_lines:
            rows = rows[rows["cell_id"].isin(list(cell_lines))]
        if time_h and "pert_itime" in rows:
            want = {f"{int(t)} h" for t in time_h}
            rows = rows[rows["pert_itime"].isin(want)]
        if min_tas is not None and "tas" in rows:
            rows = rows[rows["tas"].astype(float) >= min_tas]
        if len(rows):
            out[p] = rows
    return out


def ingest_lincs(
    gctx: str | Path,
    sig_info: str | Path,
    gene_info: str | Path,
    perts: Iterable[str],
    sig_metrics: str | Path | None = None,
    pert_info: str | Path | None = None,
    targets: dict[str, list[str]] | None = None,
    cell_lines: Iterable[str] | None = None,
    time_h: Iterable[int] | None = None,
    min_tas: float | None = None,
    gene_space: str = "bing",
    accession: str | None = None,
    n_random: int = 0,
    seed: int = 0,
) -> tuple[list[Signature], list[str]]:
    """Collapse Level 5 signatures to one consensus vector per compound.

    Consensus is the per-gene median across the selected signatures (all core
    cell lines unless ``cell_lines`` restricts it, e.g. to a patient-matched
    lineage). ``n_random`` adds that many further compounds drawn
    deterministically (``seed``) from the trt_cp pool, e.g. as a background
    for ranking. Returns (signatures, names_not_found).
    """
    si = pd.read_csv(sig_info, sep="\t", dtype=str)
    if sig_metrics:
        sm = pd.read_csv(sig_metrics, sep="\t", dtype=str)
        keep = [c for c in ("sig_id", "tas", "distil_cc_q75", "pct_self_rank_q25", "is_exemplar") if c in sm]
        si = si.merge(sm[keep], on="sig_id", how="left", suffixes=("", "_m"))
    perts = list(dict.fromkeys(p.lower() for p in perts))
    if n_random:
        pool = si[si["pert_type"] == "trt_cp"] if "pert_type" in si else si
        pool = sorted(set(pool["pert_iname"].str.lower()) - set(perts))
        rng = np.random.default_rng(seed)
        perts += [pool[i] for i in sorted(rng.choice(len(pool), min(n_random, len(pool)), replace=False))]
    groups = select_signatures(si, perts, cell_lines, time_h, min_tas)
    not_found = [p for p in perts if p not in groups]
    if not groups:
        return [], not_found

    all_ids = [s for rows in groups.values() for s in rows["sig_id"]]
    mat = read_gctx_columns(gctx, all_ids)
    gmap = load_gene_map(gene_info, gene_space)
    mat = mat.loc[mat.index.intersection(list(gmap))]
    symbols = [gmap[i] for i in mat.index]

    smiles = {}
    if pert_info:
        pi = pd.read_csv(pert_info, sep="\t", dtype=str)
        if "canonical_smiles" in pi:
            smiles = dict(zip(pi["pert_id"], pi["canonical_smiles"]))
    accession = accession or Path(gctx).name.split("_")[0]

    sigs = []
    for name, rows in groups.items():
        cons = mat[rows["sig_id"].tolist()].median(axis=1).to_numpy()
        genes, vals = dedupe_vector(symbols, cons)
        cells = sorted(rows["cell_id"].unique())
        pert_id = rows["pert_id"].mode().iat[0]
        tas = pd.to_numeric(rows.get("tas"), errors="coerce") if "tas" in rows else pd.Series(dtype=float)
        sigs.append(
            Signature(
                sig_id=name.lower(),
                kind=DRUG,
                genes=genes,
                z=vals,
                provenance=PROV_LINCS,
                modality="compound",
                tissue=cells[0] if len(cells) == 1 else "cell_line_consensus",
                sample_size=len(rows),
                contrast="LINCS L1000 Level 5 MODZ, per-gene median across signatures",
                source_accessions=[accession],
                meta={
                    "pert_id": pert_id,
                    "pert_iname": rows["pert_iname"].iat[0],
                    "cell_lines": cells,
                    "doses": sorted(rows["pert_idose"].dropna().unique().tolist()) if "pert_idose" in rows else [],
                    "time_h": sorted(rows["pert_itime"].dropna().unique().tolist()) if "pert_itime" in rows else [],
                    "n_signatures": int(len(rows)),
                    "tas_median": float(np.nanmedian(tas)) if len(tas.dropna()) else None,
                    "smiles": smiles.get(pert_id),
                    "targets": (targets or {}).get(name.lower(), []),
                    "gene_space": gene_space,
                },
            )
        )
    return sigs, not_found
