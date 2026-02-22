from __future__ import annotations

import argparse
import os
import time
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import yaml


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--in-parquet", required=True)
    parser.add_argument("--out-parquet", required=True)
    args = parser.parse_args()

    cfg = yaml.safe_load(open(args.config, "r", encoding="utf-8"))
    threshold = float(cfg["analysis"]["p_threshold_eqtl"])
    progress_every = int(cfg.get("runtime", {}).get("progress_every_row_groups", 10))

    pf = pq.ParquetFile(args.in_parquet)
    out_path = Path(args.out_parquet)
    tmp_path = Path(f"{args.out_parquet}.tmp")
    if tmp_path.exists():
        tmp_path.unlink()

    writer = None
    n_output = 0
    started = time.time()

    for rg in range(pf.metadata.num_row_groups):
        chunk = pf.read_row_group(rg).to_pandas()
        out = chunk.loc[pd.to_numeric(chunk["pvalue"], errors="coerce") < threshold].copy()
        if not out.empty:
            table = pa.Table.from_pandas(out, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(str(tmp_path), table.schema, compression="zstd")
            writer.write_table(table)
            n_output += len(out)
        if (rg + 1) % progress_every == 0:
            elapsed = time.time() - started
            print(
                f"[filter_eqtl_stream] row_groups={rg + 1}/{pf.metadata.num_row_groups} n_output={n_output} elapsed={elapsed:.1f}s",
                flush=True,
            )

    if writer is not None:
        writer.close()
    else:
        pd.DataFrame().to_parquet(str(tmp_path), index=False)

    if out_path.exists():
        out_path.unlink()
    os.replace(str(tmp_path), str(out_path))


if __name__ == "__main__":
    main()
