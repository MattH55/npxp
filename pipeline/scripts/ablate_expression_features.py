"""Does cell-line gene expression carry usable synergy information? The SynVerse test.

SynVerse (Brief Bioinform 2025) evaluated 16 published drug-synergy models and found
that none beat a naive one-hot-encoding baseline, and that models given *shuffled*
drug and cell-line features performed about as well as models given the real ones --
so their accuracy came from agent identity seen in training, not from biology. That is
the single most important result for this project, and it is reproducible here:
DrugComb supplies 740k measured drug pairs and DepMap supplies the expression.

The question this answers is narrow and practical: **if you brute-force gene
expression into a synergy model, does it buy anything?** Four feature sets, so the
comparison isolates it:

  drugs_only        drug_row and drug_col as categories. The floor: no cell context.
  identity          the above plus cell_line as a category. Cell context as pure
                    identity, learnable only for cell lines seen in training.
  expression        the drugs plus PCs of that cell line's DepMap expression, and NO
                    cell-line category. Expression must carry the cell context alone.
  expression_shuffled   the same PCs, but assigned to the WRONG cell lines. This is
                    SynVerse's ablation. Any gap between this and `expression` is what
                    real biology is worth; no gap means the features are decoration.

Three splits, because the split is what decides whether identity can be memorised:

  random            rows split at random. Every drug and cell line appears in
                    training, so identity is maximally exploitable.
  leave_cell_out    GroupKFold on cell line. Test cell lines are unseen, so their
                    identity is worthless and expression is the only route to cell
                    context. **This is the decisive split.**
  leave_drug_out    GroupKFold on drug_row. Test drugs are unseen.

Read `expression` against `expression_shuffled` under `leave_cell_out`. If they tie,
expression carries no usable information about synergy in this data, and no amount of
brute force over more expression datasets changes that.

    python scripts/ablate_expression_features.py --drugcomb ".../summary_v_1_5.csv"
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import GroupKFold, KFold

FEATURE_SETS = ["drugs_only", "identity", "identity_pair", "expression",
                "expression_shuffled", "identity_pair_plus_expression"]
SPLITS = ["random", "leave_cell_out", "leave_drug_out"]


def load_expression(depmap: Path, n_pcs: int = 50) -> tuple[pd.DataFrame, dict[str, str]]:
    """Cell line x PCs of DepMap expression, and a name -> model-id map."""
    model = pd.read_csv(depmap / "Model.csv", low_memory=False)
    expr = pd.read_csv(depmap / "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                       index_col=0, low_memory=False)
    expr = expr.dropna(axis=1, how="any")
    # most variable genes first, so the PCA is not dominated by noise floor
    keep = expr.var().nlargest(4000).index
    P = PCA(n_components=min(n_pcs, min(expr.shape) - 1), random_state=0)
    pcs = pd.DataFrame(P.fit_transform(expr[keep].to_numpy()), index=expr.index,
                       columns=[f"pc{i}" for i in range(P.n_components_)])
    name_to_id = {}
    for col in ("StrippedCellLineName", "CellLineName", "ModelID"):
        if col in model.columns:
            for nm, mid in zip(model[col].astype(str), model["ModelID"].astype(str)):
                key = "".join(ch for ch in nm.upper() if ch.isalnum())
                if key and key not in name_to_id:
                    name_to_id[key] = mid
    return pcs, name_to_id


def score(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    ok = np.isfinite(y_true) & np.isfinite(y_pred)
    if ok.sum() < 50 or np.std(y_pred[ok]) == 0:
        return {"pearson": float("nan"), "spearman": float("nan")}
    a, b = y_true[ok], y_pred[ok]
    pear = float(np.corrcoef(a, b)[0, 1])
    ra = pd.Series(a).rank().to_numpy()
    rb = pd.Series(b).rank().to_numpy()
    return {"pearson": pear, "spearman": float(np.corrcoef(ra, rb)[0, 1])}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--drugcomb", required=True)
    ap.add_argument("--depmap", default="data/raw/depmap")
    ap.add_argument("--target", default="synergy_zip",
                    choices=["synergy_zip", "synergy_bliss", "synergy_loewe", "synergy_hsa"])
    ap.add_argument("--max-rows", type=int, default=150_000)
    ap.add_argument("--n-pcs", type=int, default=50)
    ap.add_argument("--folds", type=int, default=4)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="out/expression_ablation")
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    d = pd.read_csv(a.drugcomb, low_memory=False,
                    usecols=["drug_row", "drug_col", "cell_line_name", a.target])
    d = d.dropna(subset=["drug_row", "drug_col", "cell_line_name", a.target])
    print(f"DrugComb: {len(d):,} rows with {a.target}", file=sys.stderr)

    pcs, name_to_id = load_expression(Path(a.depmap), a.n_pcs)
    d["_key"] = ["".join(c for c in str(n).upper() if c.isalnum()) for n in d["cell_line_name"]]
    d["_model"] = d["_key"].map(name_to_id)
    d = d[d["_model"].isin(pcs.index)]
    print(f"{len(d):,} rows on {d['_model'].nunique()} cell lines matched to DepMap "
          f"expression", file=sys.stderr)
    if len(d) > a.max_rows:
        d = d.sample(a.max_rows, random_state=a.seed)
    d = d.reset_index(drop=True)

    # Identity is encoded as a fold-safe target encoding rather than a category, because
    # there are ~1,445 distinct drugs and sklearn's histogram GBT caps a categorical at
    # 255 levels. A per-level mean of the target IS the identity channel -- it is what a
    # one-hot tree would learn about that level -- and it is computed inside each fold
    # from the training rows only, so no test information reaches the encoding.
    ident_cols = ["drug_row", "drug_col", "cell_line_name"]
    keys = d[ident_cols].astype(str)
    # Pair and triple keys. A per-level mean captures only MAIN effects, and the
    # accuracy published synergy models report comes largely from memorising the PAIR
    # (and the pair in that cell line). Without these the identity channel is too weak
    # to stand in for SynVerse's one-hot baseline.
    keys["pair"] = np.where(keys["drug_row"] < keys["drug_col"],
                            keys["drug_row"] + "|" + keys["drug_col"],
                            keys["drug_col"] + "|" + keys["drug_row"])
    keys["triple"] = keys["pair"] + "|" + keys["cell_line_name"]
    X_pcs = pcs.loc[d["_model"]].to_numpy()
    # Shuffle which cell line's expression each row receives, by permuting the MAP
    # from cell line to profile. Every profile stays a real profile; only the pairing
    # with the row's actual cell line is destroyed.
    cells = d["_model"].unique()
    permuted = dict(zip(cells, rng.permutation(cells)))
    X_pcs_shuf = pcs.loc[d["_model"].map(permuted)].to_numpy()
    y = d[a.target].to_numpy(float)

    def encode(train_idx: np.ndarray, cols: list[str]) -> np.ndarray:
        """Per-level mean of the target, learned on `train_idx` only."""
        out = np.empty((len(d), len(cols)))
        for j, c in enumerate(cols):
            grand = float(y[train_idx].mean())
            means = pd.Series(y[train_idx]).groupby(keys[c].to_numpy()[train_idx]).mean()
            out[:, j] = keys[c].map(means).fillna(grand).to_numpy(float)
        return out

    # (identity columns used, whether to append expression, which expression matrix)
    designs = {
        "drugs_only": (["drug_row", "drug_col"], None),
        "identity": (ident_cols, None),
        "identity_pair": (ident_cols + ["pair", "triple"], None),
        "expression": (["drug_row", "drug_col"], X_pcs),
        "expression_shuffled": (["drug_row", "drug_col"], X_pcs_shuf),
        "identity_pair_plus_expression": (ident_cols + ["pair", "triple"], X_pcs),
    }
    groups = {"random": None, "leave_cell_out": d["_model"].to_numpy(),
              "leave_drug_out": d["drug_row"].astype(str).to_numpy()}

    rows = []
    for split, g in groups.items():
        splitter = (KFold(a.folds, shuffle=True, random_state=a.seed) if g is None
                    else GroupKFold(a.folds))
        folds = list(splitter.split(d, groups=g) if g is not None else splitter.split(d))
        for name, (cols, extra) in designs.items():
            preds = np.full(len(d), np.nan)
            for tr, te in folds:
                E = encode(tr, cols)
                X = E if extra is None else np.hstack([E, extra])
                m = HistGradientBoostingRegressor(
                    max_iter=200, learning_rate=0.1, random_state=a.seed)
                m.fit(X[tr], y[tr])
                preds[te] = m.predict(X[te])
            s = score(y, preds)
            rows.append({"split": split, "features": name, **s})
            print(f"  {split:15s} {name:22s} pearson {s['pearson']:+.3f}  "
                  f"spearman {s['spearman']:+.3f}", file=sys.stderr)

    t = pd.DataFrame(rows)
    t.to_csv(out / "ablation.tsv", sep="\t", index=False)
    piv = t.pivot(index="split", columns="features", values="pearson")[FEATURE_SETS]
    print("\nPearson r by feature set and split:\n")
    print(piv.to_string(float_format=lambda v: f"{v:+.3f}"))

    verdict = {
        "target": a.target, "n_rows": int(len(d)),
        "n_cell_lines": int(d["_model"].nunique()),
        "pearson": piv.to_dict(),
        "expression_vs_shuffled": {
            s: round(float(piv.loc[s, "expression"] - piv.loc[s, "expression_shuffled"]), 4)
            for s in piv.index},
        "expression_vs_identity": {
            s: round(float(piv.loc[s, "expression"] - piv.loc[s, "identity"]), 4)
            for s in piv.index},
        "reading":
            "Compare expression against expression_shuffled under leave_cell_out. That "
            "split makes a test cell line's identity worthless, so expression is the "
            "only route to cell context; if the two tie, expression carries no usable "
            "synergy information and brute-forcing more expression data cannot help.",
    }
    (out / "verdict.json").write_text(json.dumps(verdict, indent=2, default=float))
    gap = verdict["expression_vs_shuffled"].get("leave_cell_out")
    print(f"\nreal minus shuffled expression, leave_cell_out: {gap:+.4f}")
    print("Expression carries usable information." if gap is not None and gap > 0.02
          else "Expression adds nothing over shuffled features: it is decoration here.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
