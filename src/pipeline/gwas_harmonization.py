from __future__ import annotations

import argparse
import os
import time
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import yaml

from src.utils.genetics import normalize_chr


def table_params(path: Path) -> tuple[str, str | None]:
    name = path.name.lower()
    if name.endswith(".csv") or name.endswith(".csv.gz"):
        sep = ","
    else:
        sep = "\t"
    compression = "gzip" if name.endswith(".gz") else None
    return sep, compression


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--out-parquet", required=True)
    parser.add_argument("--qc", required=True)
    args = parser.parse_args()

    cfg = yaml.safe_load(open(args.config, "r", encoding="utf-8"))
    manifest = yaml.safe_load(open(args.manifest, "r", encoding="utf-8"))

    target_build = cfg.get("genome_build", "GRCh38")
    chunksize = int(cfg.get("runtime", {}).get("chunksize", 250000))
    writer = None
    qc_rows: list[dict] = []
    out_path = Path(args.out_parquet)
    tmp_path = Path(f"{args.out_parquet}.tmp")
    if tmp_path.exists():
        tmp_path.unlink()

    for trait in manifest["traits"]:
        path = Path(trait["file_path"])
        mandatory = bool(trait.get("is_mandatory", True))
        if not path.exists():
            qc_rows.append(
                {
                    "trait_id": trait["trait_id"],
                    "status": "missing_mandatory" if mandatory else "missing_optional",
                    "path": str(path),
                    "n_rows": 0,
                    "build": trait.get("build", ""),
                }
            )
            continue

        if trait.get("build", target_build) != target_build:
            qc_rows.append(
                {
                    "trait_id": trait["trait_id"],
                    "status": "build_mismatch_no_liftover",
                    "path": str(path),
                    "n_rows": 0,
                    "build": trait.get("build", ""),
                }
            )
            continue

        sep, compression = table_params(path)
        cols = [
            trait["chr_col"],
            trait["pos_col"],
            trait["ea_col"],
            trait["oa_col"],
            trait["beta_col"],
            trait["se_col"],
            trait["p_col"],
            trait["rsid_col"],
        ]
        n_rows = 0
        n_chunks = 0
        started = time.time()
        progress_every = int(cfg.get("runtime", {}).get("progress_every_chunks", 25))
        print(f"[gwas_harmonization] start trait={trait['trait_id']} path={path}", flush=True)
        schema_ok = True
        try:
            reader = pd.read_csv(
                path,
                sep=sep,
                compression=compression,
                usecols=cols,
                chunksize=chunksize,
                low_memory=False,
            )
            for raw in reader:
                n_chunks += 1
                frame = pd.DataFrame(
                    {
                        "trait_id": trait["trait_id"],
                        "trait_label": trait["trait_label"],
                        "ancestry": trait.get("ancestry", "NA"),
                        "chr": raw[trait["chr_col"]].astype(str).map(normalize_chr),
                        "pos": pd.to_numeric(raw[trait["pos_col"]], errors="coerce"),
                        "effect_allele": raw[trait["ea_col"]].astype(str).str.upper(),
                        "other_allele": raw[trait["oa_col"]].astype(str).str.upper(),
                        "beta": pd.to_numeric(raw[trait["beta_col"]], errors="coerce"),
                        "se": pd.to_numeric(raw[trait["se_col"]], errors="coerce"),
                        "p": pd.to_numeric(raw[trait["p_col"]], errors="coerce"),
                        "rsid": raw[trait["rsid_col"]].astype(str),
                    }
                )
                frame = frame.dropna(subset=["pos", "beta", "se", "p"]).copy()
                frame = frame.dropna(subset=["chr"]).copy()
                if frame.empty:
                    continue
                table = pa.Table.from_pandas(frame, preserve_index=False)
                if writer is None:
                    writer = pq.ParquetWriter(str(tmp_path), table.schema)
                writer.write_table(table)
                n_rows += len(frame)
                if n_chunks % progress_every == 0:
                    elapsed = time.time() - started
                    rate = n_rows / elapsed if elapsed > 0 else 0.0
                    print(
                        f"[gwas_harmonization] trait={trait['trait_id']} chunks={n_chunks} rows={n_rows} rate_rows_s={rate:.1f}",
                        flush=True,
                    )
        except ValueError:
            schema_ok = False

        if not schema_ok:
            qc_rows.append(
                {
                    "trait_id": trait["trait_id"],
                    "status": "schema_mismatch",
                    "path": str(path),
                    "n_rows": 0,
                    "build": trait.get("build", ""),
                }
            )
            continue

        qc_rows.append(
            {
                "trait_id": trait["trait_id"],
                "status": "ok",
                "path": str(path),
                "n_rows": n_rows,
                "build": trait.get("build", ""),
            }
        )
        elapsed = time.time() - started
        print(
            f"[gwas_harmonization] done trait={trait['trait_id']} chunks={n_chunks} rows={n_rows} elapsed={elapsed:.1f}s",
            flush=True,
        )

    if writer is not None:
        writer.close()
    else:
        out = pd.DataFrame(
            columns=[
                "trait_id",
                "trait_label",
                "ancestry",
                "chr",
                "pos",
                "effect_allele",
                "other_allele",
                "beta",
                "se",
                "p",
                "rsid",
            ]
        )
        out.to_parquet(str(tmp_path), index=False)
    if out_path.exists():
        out_path.unlink()
    os.replace(str(tmp_path), str(out_path))
    pd.DataFrame(qc_rows).to_csv(args.qc, sep="\t", index=False)


if __name__ == "__main__":
    main()
