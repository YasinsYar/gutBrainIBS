from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Iterable

import pandas as pd


def ensure_parent(path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def ensure_dir(path: str | Path) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def sha256_file(path: str | Path, chunk_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(chunk_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def read_tsv_head(path: str | Path, n: int = 5) -> list[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    lines: list[str] = []
    with opener(path, "rt", encoding="utf-8", errors="ignore") as handle:
        for _ in range(n):
            line = handle.readline()
            if not line:
                break
            lines.append(line.rstrip("\n"))
    return lines


def write_json(data: dict, path: str | Path) -> None:
    ensure_parent(path)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)


def write_tsv(df: pd.DataFrame, path: str | Path) -> None:
    ensure_parent(path)
    df.to_csv(path, sep="\t", index=False)


def append_manifest_rows(path: str | Path, rows: Iterable[dict]) -> None:
    ensure_parent(path)
    rows = list(rows)
    if not rows:
        return
    frame = pd.DataFrame(rows)
    if Path(path).exists():
        frame.to_csv(path, sep="\t", index=False, mode="a", header=False)
    else:
        frame.to_csv(path, sep="\t", index=False)
