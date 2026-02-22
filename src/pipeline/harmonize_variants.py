from __future__ import annotations

import argparse
import os
import time
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import yaml

from src.utils.genetics import is_palindromic


def process_chunk(frame: pd.DataFrame, drop_palindromic: bool) -> tuple[pd.DataFrame, int, int]:
    if frame.empty:
        return frame, 0, 0

    frame = frame.copy()
    frame["ref"] = frame["ref"].astype(str).str.upper()
    frame["alt"] = frame["alt"].astype(str).str.upper()
    frame["is_palindromic"] = [is_palindromic(r, a) for r, a in zip(frame["ref"], frame["alt"])]
    frame["effect_allele"] = frame["alt"]
    frame["other_allele"] = frame["ref"]
    frame["was_flipped"] = False

    n_pal = int(frame["is_palindromic"].sum())
    if drop_palindromic:
        frame = frame.loc[~frame["is_palindromic"]].copy()

    before = len(frame)
    frame = frame.drop_duplicates(subset=["chr", "pos", "ref", "alt", "gene_id", "dataset_id", "tissue"])
    n_dups = before - len(frame)
    return frame, n_pal, n_dups


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--in-parquet", required=True)
    parser.add_argument("--out-parquet", required=True)
    parser.add_argument("--qc", required=True)
    parser.add_argument("--exclude-palindromic", choices=["true", "false"], default=None)
    args = parser.parse_args()

    cfg = yaml.safe_load(open(args.config, "r", encoding="utf-8"))
    if args.exclude_palindromic is None:
        drop_pal = bool(cfg["analysis"].get("exclude_palindromic", True))
    else:
        drop_pal = args.exclude_palindromic.lower() == "true"
    progress_every = int(cfg.get("runtime", {}).get("progress_every_row_groups", 10))

    pf = pq.ParquetFile(args.in_parquet)
    n_input = int(pf.metadata.num_rows)

    out_path = Path(args.out_parquet)
    tmp_path = Path(f"{args.out_parquet}.tmp")
    if tmp_path.exists():
        tmp_path.unlink()

    writer = None
    n_output = 0
    n_pal = 0
    n_dups = 0
    pre_dedup_total = 0
    started = time.time()

    for rg in range(pf.metadata.num_row_groups):
        chunk = pf.read_row_group(rg).to_pandas()
        pre_dedup_total += len(chunk)
        out, pal_count, dup_count = process_chunk(chunk, drop_pal)
        n_pal += pal_count
        n_dups += dup_count
        if not out.empty:
            table = pa.Table.from_pandas(out, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(str(tmp_path), table.schema, compression="zstd")
            writer.write_table(table)
            n_output += len(out)
        if (rg + 1) % progress_every == 0:
            elapsed = time.time() - started
            print(
                f"[harmonize_stream] row_groups={rg + 1}/{pf.metadata.num_row_groups} n_output={n_output} elapsed={elapsed:.1f}s",
                flush=True,
            )

    if writer is not None:
        writer.close()
    else:
        pd.DataFrame().to_parquet(str(tmp_path), index=False)

    if out_path.exists():
        out_path.unlink()
    os.replace(str(tmp_path), str(out_path))

    qc = pd.DataFrame(
        [
            {"metric": "n_input", "value": n_input},
            {"metric": "n_palindromic", "value": n_pal},
            {"metric": "n_deduplicated", "value": n_dups},
            {"metric": "n_output", "value": n_output},
            {"metric": "palindrome_drop_rate", "value": (n_pal / n_input) if n_input else 0.0},
            {"metric": "duplicate_drop_rate", "value": (n_dups / max(pre_dedup_total, 1))},
        ]
    )
    qc.to_csv(args.qc, sep="\t", index=False)


if __name__ == "__main__":
    main()
