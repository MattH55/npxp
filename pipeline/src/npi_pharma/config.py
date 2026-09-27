"""Scoring configuration loading (bundled defaults + user overrides)."""

from __future__ import annotations

import copy
from importlib import resources
from pathlib import Path
from typing import Any

import yaml


def _merge(base: dict, over: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge(out[k], v)
        else:
            out[k] = v
    return out


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    base = yaml.safe_load(resources.files("npi_pharma").joinpath("data", "scoring.yaml").read_text())
    if path:
        base = _merge(base, yaml.safe_load(Path(path).read_text()) or {})
    return base


def canonical_tissue(tissue: str | None, cfg: dict[str, Any]) -> str | None:
    if tissue is None:
        return None
    t = str(tissue).strip().lower().replace(" ", "_").replace("-", "_")
    for canon, syns in cfg.get("tissues", {}).items():
        if t == canon or t in syns:
            return canon
    return t
