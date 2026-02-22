from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import time
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
import yaml

from src.utils.duckdb_utils import connect_duckdb, sql_path, sql_quote
from src.utils.io import ensure_parent


BASE_COLUMNS = ["chromosome", "position", "ref", "alt", "variant", "pvalue", "beta", "se"]
GENE_CANDIDATES = ["gene_id", "molecular_trait_object_id", "molecular_trait_id"]
OPTIONAL_COLUMNS = ["maf", "rsid"]


def iter_input_files(colon_manifest: Path, brain_manifest: Path) -> list[dict]:
    files: list[dict] = []
    for manifest in [colon_manifest, brain_manifest]:
        if not manifest.exists():
            continue
        frame = pd.read_csv(manifest, sep="\t")
        for row in frame.itertuples(index=False):
            status = str(getattr(row, "status", ""))
            if status.startswith("failed"):
                continue
            path = Path(getattr(row, "local_path"))
            if not path.exists():
                continue
            files.append(
                {
                    "dataset_id": str(getattr(row, "dataset_id")),
                    "tissue": str(getattr(row, "tissue")),
                    "source": str(getattr(row, "source", "eqtl_catalogue")),
                    "path": path,
                }
            )
    return files


def detect_columns(path: Path) -> set[str]:
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as handle:
        header = handle.readline().strip().split("\t")
    return set(header)


def part_filename(item: dict) -> str:
    key = f"{item['dataset_id']}|{item['tissue']}|{Path(item['path']).resolve()}"
    digest = hashlib.sha1(key.encode("utf-8")).hexdigest()[:12]
    return f"{item['dataset_id']}_{digest}.parquet"


def part_is_fresh(part_path: Path, source_path: Path) -> bool:
    if not part_path.exists() or not source_path.exists():
        return False
    try:
        return part_path.stat().st_mtime >= source_path.stat().st_mtime
    except OSError:
        return False


def parse_dataset_to_part(con, item: dict, part_path: Path, row_group_size: int) -> tuple[int, int, str]:
    cols = detect_columns(item["path"])
    missing = [c for c in BASE_COLUMNS if c not in cols]
    if missing:
        return 0, 0, f"missing_columns:{','.join(missing)}"

    gene_col = next((c for c in GENE_CANDIDATES if c in cols), None)

    gene_expr = f"regexp_replace(CAST({gene_col} AS VARCHAR), '\\\\..*$', '')" if gene_col else "''"
    maf_expr = "TRY_CAST(maf AS DOUBLE)" if "maf" in cols else "NULL::DOUBLE"
    rsid_expr = "COALESCE(CAST(rsid AS VARCHAR), '')" if "rsid" in cols else "''"

    tmp_part = Path(f"{part_path}.tmp")
    if tmp_part.exists():
        tmp_part.unlink()

    started = time.time()
    print(f"[parse_eqtl_duckdb] start dataset={item['dataset_id']} tissue={item['tissue']}", flush=True)
    sql = f"""
    COPY (
      SELECT
        {sql_quote(item['source'])}::VARCHAR AS source,
        {sql_quote(item['dataset_id'])}::VARCHAR AS dataset_id,
        {sql_quote(item['tissue'])}::VARCHAR AS tissue,
        CASE WHEN starts_with(CAST(chromosome AS VARCHAR), 'chr')
             THEN CAST(chromosome AS VARCHAR)
             ELSE 'chr' || CAST(chromosome AS VARCHAR) END AS chr,
        TRY_CAST(position AS BIGINT) AS pos,
        upper(CAST(ref AS VARCHAR)) AS ref,
        upper(CAST(alt AS VARCHAR)) AS alt,
        CAST(variant AS VARCHAR) AS variant,
        {gene_expr} AS gene_id,
        TRY_CAST(beta AS DOUBLE) AS beta,
        TRY_CAST(se AS DOUBLE) AS se,
        TRY_CAST(pvalue AS DOUBLE) AS pvalue,
        {maf_expr} AS maf,
        {rsid_expr} AS rsid
      FROM read_csv_auto(
        {sql_path(item['path'])},
        delim='\\t',
        header=TRUE,
        compression='gzip',
        all_varchar=TRUE,
        ignore_errors=TRUE
      )
      WHERE TRY_CAST(position AS BIGINT) IS NOT NULL
        AND TRY_CAST(beta AS DOUBLE) IS NOT NULL
        AND TRY_CAST(se AS DOUBLE) IS NOT NULL
        AND TRY_CAST(pvalue AS DOUBLE) IS NOT NULL
    ) TO {sql_path(tmp_part)} (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE {int(row_group_size)});
    """
    con.execute(sql)

    if part_path.exists():
        part_path.unlink()
    os.replace(str(tmp_part), str(part_path))
    pf = pq.ParquetFile(part_path)
    elapsed = time.time() - started
    n_rows = int(pf.metadata.num_rows)
    n_groups = int(pf.metadata.num_row_groups)
    print(
        f"[parse_eqtl_duckdb] done dataset={item['dataset_id']} rows={n_rows} row_groups={n_groups} elapsed={elapsed:.1f}s",
        flush=True,
    )
    return n_rows, n_groups, "parsed"


def clear_stale_tmp(temp_dir: Path) -> None:
    if not temp_dir.exists():
        return
    for p in temp_dir.glob("duckdb_temp_storage_*.tmp"):
        try:
            p.unlink()
        except OSError:
            continue


def combine_parts(part_paths: list[Path], out_path: Path, delete_parts: bool) -> int:
    existing = [p for p in part_paths if p.exists()]
    total_rows = 0
    for p in existing:
        total_rows += int(pq.ParquetFile(p).metadata.num_rows)

    tmp_out = Path(f"{out_path}.tmp")
    if tmp_out.exists():
        tmp_out.unlink()

    if not existing:
        pd.DataFrame(
            columns=[
                "source",
                "dataset_id",
                "tissue",
                "chr",
                "pos",
                "ref",
                "alt",
                "variant",
                "gene_id",
                "beta",
                "se",
                "pvalue",
                "maf",
                "rsid",
            ]
        ).to_parquet(str(tmp_out), index=False)
    else:
        writer = None
        copied_rows = 0
        for idx, part in enumerate(existing, start=1):
            pf = pq.ParquetFile(part)
            for rg in range(pf.metadata.num_row_groups):
                table = pf.read_row_group(rg)
                if writer is None:
                    writer = pq.ParquetWriter(str(tmp_out), table.schema, compression="zstd")
                writer.write_table(table)
                copied_rows += table.num_rows
            print(
                f"[parse_eqtl_duckdb] combine part={idx}/{len(existing)} rows_copied={copied_rows}",
                flush=True,
            )
            if delete_parts:
                try:
                    part.unlink()
                    print(f"[parse_eqtl_duckdb] removed_part={part}", flush=True)
                except OSError:
                    pass
        if writer is not None:
            writer.close()

    if out_path.exists():
        out_path.unlink()
    os.replace(str(tmp_out), str(out_path))
    return total_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--colon-manifest", required=True)
    parser.add_argument("--brain-manifest", required=True)
    parser.add_argument("--out-parquet", required=True)
    parser.add_argument("--out-manifest", required=True)
    parser.add_argument("--force-reparse", action="store_true")
    args = parser.parse_args()

    cfg = yaml.safe_load(open(args.config, "r", encoding="utf-8"))
    files = iter_input_files(Path(args.colon_manifest), Path(args.brain_manifest))

    ensure_parent(args.out_parquet)
    out_path = Path(args.out_parquet)
    parts_dir = Path("data_interim/eqtl_parts")
    parts_dir.mkdir(parents=True, exist_ok=True)
    parse_rows: list[dict] = []
    part_paths: list[Path] = []
    row_group_size = int(cfg.get("runtime", {}).get("parquet_row_group_size", 500000))
    delete_parts = bool(cfg.get("runtime", {}).get("cleanup_eqtl_parts_after_combine", True))
    duckdb_tmp = Path(cfg.get("runtime", {}).get("duckdb_temp_dir", "data_interim/duckdb_tmp"))
    clear_stale_tmp(duckdb_tmp)

    con = connect_duckdb(cfg)
    try:
        for item in files:
            part_path = parts_dir / part_filename(item)
            part_paths.append(part_path)
            if part_is_fresh(part_path, Path(item["path"])) and not args.force_reparse:
                pf = pq.ParquetFile(part_path)
                parse_rows.append(
                    {
                        "dataset_id": item["dataset_id"],
                        "tissue": item["tissue"],
                        "path": str(item["path"]),
                        "rows_parsed": int(pf.metadata.num_rows),
                        "chunks_parsed": int(pf.metadata.num_row_groups),
                        "status": "skipped_existing",
                        "part_path": str(part_path),
                    }
                )
                print(f"[parse_eqtl_duckdb] skip dataset={item['dataset_id']} existing_part={part_path}", flush=True)
                continue

            rows, chunks, status = parse_dataset_to_part(con, item, part_path, row_group_size)
            parse_rows.append(
                {
                    "dataset_id": item["dataset_id"],
                    "tissue": item["tissue"],
                    "path": str(item["path"]),
                    "rows_parsed": rows,
                    "chunks_parsed": chunks,
                    "status": status,
                    "part_path": str(part_path),
                }
            )

        total_rows = combine_parts(part_paths, out_path, delete_parts=delete_parts)
    finally:
        con.close()

    print(f"[parse_eqtl_duckdb] combined_parts={len(part_paths)} total_rows={total_rows}", flush=True)
    pd.DataFrame(parse_rows).to_csv(args.out_manifest, sep="\t", index=False)


if __name__ == "__main__":
    main()
