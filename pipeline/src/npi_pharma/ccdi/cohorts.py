"""Pediatric tumour / model cohorts from OpenPedCan, and model-to-tumour concordance."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

TUMOUR, MODEL = "tumour", "model"
MODEL_COMPOSITIONS = {"Derived Cell Line", "Patient Derived Xenograft"}


def load_diagnosis_map(path: str | Path) -> pd.DataFrame:
    m = pd.read_csv(path, sep="\t", comment="#", dtype=str).fillna("")
    return m.set_index("cancer_group")


def map_groups(groups: list[str], dmap: pd.DataFrame) -> tuple[dict[str, str], list[str]]:
    """cancer_group -> OncoTree code; returns (mapped, unmapped incl. intentionally blank)."""
    mapped, unmapped = {}, []
    for g in sorted(set(groups)):
        code = dmap["oncotree_code"].get(g, "")
        (mapped.__setitem__(g, code) if code else unmapped.append(g))
    return mapped, unmapped


def tpm_to_parquet(rds: str | Path, out: str | Path) -> Path:
    """Collapsed RSEM TPM .rds -> log2(TPM+1) float32 parquet (genes x samples)."""
    from .openpedcan import read_tpm_rds

    df = read_tpm_rds(rds)
    df = np.log2(df.clip(lower=0) + 1).astype(np.float32)
    df = df[~df.index.duplicated()]
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out)
    return out


def rna_samples(hist: pd.DataFrame, expr_columns, independent: pd.DataFrame | None = None) -> pd.DataFrame:
    """RNA-seq biospecimens with expression, labelled tumour or model.

    Tumours: sample_type Tumor, solid tissue / blood / marrow, restricted to the
    independent primary specimens (one per patient) when given. Models: derived
    cell lines and PDXs (all kept; they are distinct models).
    """
    h = hist[(hist["experimental_strategy"] == "RNA-Seq") & (hist["sample_type"] == "Tumor")]
    h = h[h["Kids_First_Biospecimen_ID"].isin(set(expr_columns)) & h["cancer_group"].notna()]
    h = h[h["cohort"] != "TCGA"]  # adult; kept out of pediatric cohorts
    kind = np.where(h["composition"].isin(MODEL_COMPOSITIONS), MODEL, TUMOUR)
    h = h.assign(kind=kind)
    if independent is not None:
        keep = set(independent["Kids_First_Biospecimen_ID"])
        h = h[(h["kind"] == MODEL) | h["Kids_First_Biospecimen_ID"].isin(keep)]
    return h


def cohort_counts(samples: pd.DataFrame) -> pd.DataFrame:
    t = samples.pivot_table(index="cancer_group", columns="kind", values="Kids_First_Biospecimen_ID",
                            aggfunc="count", fill_value=0)
    for k in (TUMOUR, MODEL):
        if k not in t:
            t[k] = 0
    return t[[TUMOUR, MODEL]].sort_values(TUMOUR, ascending=False)


def _spearman_rows(a: np.ndarray, centroid: np.ndarray) -> np.ndarray:
    ra = np.argsort(np.argsort(a, axis=1), axis=1).astype(float)
    rc = np.argsort(np.argsort(centroid)).astype(float)
    ra -= ra.mean(axis=1, keepdims=True)
    rc -= rc.mean()
    return (ra @ rc) / (np.linalg.norm(ra, axis=1) * np.linalg.norm(rc))


def concordance(
    expr: pd.DataFrame,
    samples: pd.DataFrame,
    min_tumours: int = 20,
    n_var: int = 1000,
    min_tpm_log: float = 1.0,
) -> pd.DataFrame:
    """Model -> tumour-group concordance, with a cross-group null.

    Run separately per model type (cell line, PDX). For each type the reference
    groups are the cancer_groups with >= ``min_tumours`` tumours AND >= 1 model
    of that type, so tumours and models are compared over the same set of
    groups. Within that set:

    - genes: the ``n_var`` expressed genes with the highest variance between
      the groups' tumour centroids;
    - tumour centroids are centred on their mean, and models on the mean of the
      per-group model means. Separate centring cancels the global
      tumour-vs-model offset (stroma, immune infiltrate, culture or mouse host);
      it borrows Celligner's idea and is not Celligner;
    - concordance = Spearman(model, own-group centroid); null_percentile = the
      share of same-type models from other groups with a lower correlation to
      that centroid; best_group = the reference centroid the model matches best.

    Models whose group is not a reference group get ``no_cohort`` (best_group
    is still reported).
    """
    tum = samples[samples["kind"] == TUMOUR]
    mod = samples[(samples["kind"] == MODEL) & samples["Kids_First_Biospecimen_ID"].isin(expr.columns)]
    big = {g for g, n in tum["cancer_group"].value_counts().items() if n >= min_tumours}
    rows = []
    for comp, info in mod.groupby("composition"):
        info = info.set_index("Kids_First_Biospecimen_ID")
        ref = sorted(big & set(info["cancer_group"]))
        if len(ref) < 2:
            for m, r in info.iterrows():
                rows.append({"biospecimen_id": m, "cancer_group": r["cancer_group"], "composition": comp,
                             "n_reference_groups": len(ref), "flag": "no_cohort"})
            continue
        cent = pd.DataFrame({g: expr[tum.loc[tum["cancer_group"] == g, "Kids_First_Biospecimen_ID"]].mean(axis=1)
                             for g in ref})
        expressed = cent.index[(cent.max(axis=1) >= min_tpm_log).to_numpy()]
        genes = cent.loc[expressed].var(axis=1).sort_values(ascending=False).index[:n_var]
        C = cent.loc[genes]
        C = C.sub(C.mean(axis=1), axis=0)
        ids = info.index.tolist()
        M = expr.loc[genes, ids]
        in_ref = info["cancer_group"].isin(ref).to_numpy()
        gm = M.loc[:, in_ref].T.groupby(info.loc[in_ref, "cancer_group"].values).mean().T
        M = M.sub(gm.mean(axis=1), axis=0)
        S = pd.DataFrame({g: _spearman_rows(M.T.to_numpy(dtype=float), C[g].to_numpy(dtype=float)) for g in ref},
                         index=ids)
        for m in ids:
            g = info.at[m, "cancer_group"]
            base = {"biospecimen_id": m, "cancer_group": g, "composition": comp, "n_reference_groups": len(ref),
                    "n_tumours": int((tum["cancer_group"] == g).sum()),
                    "best_group": S.loc[m].idxmax(), "best_concordance": float(S.loc[m].max())}
            if g not in ref:
                rows.append(base | {"flag": "no_cohort"})
                continue
            own = S.at[m, g]
            null = S.loc[info.index[(info["cancer_group"] != g).to_numpy()], g]
            pct = float((null < own).mean() * 100) if len(null) else np.nan
            flag = "concordant" if (pct >= 95 and base["best_group"] == g) else (
                "discordant" if pct < 50 else "intermediate")
            rows.append(base | {"concordance": float(own), "null_percentile": pct, "flag": flag})
    out = pd.DataFrame(rows)
    meta = mod.set_index("Kids_First_Biospecimen_ID")[["sample_id", "cohort"]]
    return out.merge(meta, left_on="biospecimen_id", right_index=True, how="left")
