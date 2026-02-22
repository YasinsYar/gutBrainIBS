from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


PALINDROMIC = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}


def strip_ensg_version(gene_id: str) -> str:
    if not isinstance(gene_id, str):
        return gene_id
    return gene_id.split(".")[0]


def normalize_chr(chrom: str | int) -> str:
    if chrom is None:
        return np.nan
    if pd.isna(chrom):
        return np.nan
    value = str(chrom).strip()
    if not value or value.lower() in {"nan", "none", "na"}:
        return np.nan
    value = value.replace("chr", "")
    if value.endswith(".0") and value[:-2].isdigit():
        value = value[:-2]
    return f"chr{value}"


def is_palindromic(ref: str, alt: str) -> bool:
    return (str(ref).upper(), str(alt).upper()) in PALINDROMIC


def add_variant_key(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["variant"] = (
        out["chr"].astype(str)
        + "_"
        + out["pos"].astype(str)
        + "_"
        + out["ref"].astype(str)
        + "_"
        + out["alt"].astype(str)
    )
    return out


def clump_by_distance(df: pd.DataFrame, p_col: str, chr_col: str, pos_col: str, kb: int) -> pd.DataFrame:
    if df.empty:
        return df
    selected = []
    window = kb * 1000
    for chrom, sub in df.sort_values([chr_col, p_col]).groupby(chr_col):
        taken = []
        for row in sub.itertuples(index=False):
            pos = getattr(row, pos_col)
            if any(abs(pos - t) <= window for t in taken):
                continue
            selected.append(row)
            taken.append(pos)
    return pd.DataFrame(selected, columns=df.columns)


def safe_read_parquet(path: str | Path) -> pd.DataFrame:
    if not Path(path).exists():
        return pd.DataFrame()
    return pd.read_parquet(path)
