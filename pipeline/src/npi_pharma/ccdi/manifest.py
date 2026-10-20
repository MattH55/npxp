"""Download provenance: MANIFEST.json per source directory."""

from __future__ import annotations

import datetime as dt
import hashlib
import json
from pathlib import Path
from typing import Any


def file_digest(path: str | Path, algo: str = "sha256", chunk: int = 1 << 20) -> str:
    h = hashlib.new(algo)
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def now_utc() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load(dir_: str | Path) -> dict[str, Any]:
    p = Path(dir_) / "MANIFEST.json"
    return json.loads(p.read_text()) if p.exists() else {}


def record(dir_: str | Path, source: dict[str, Any], files: dict[str, dict[str, Any]]) -> Path:
    """Merge ``files`` (name -> {url, sha256, size, ...}) into the directory's manifest."""
    d = Path(dir_)
    d.mkdir(parents=True, exist_ok=True)
    m = load(d)
    m.update({k: v for k, v in source.items()})
    m.setdefault("files", {}).update(files)
    m["updated"] = now_utc()
    p = d / "MANIFEST.json"
    p.write_text(json.dumps(m, indent=2, sort_keys=True) + "\n")
    return p
