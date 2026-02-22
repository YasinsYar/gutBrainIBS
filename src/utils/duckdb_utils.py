from __future__ import annotations

from pathlib import Path

import duckdb


def sql_quote(value: str) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def sql_path(path: str | Path) -> str:
    return sql_quote(str(Path(path).resolve().as_posix()))


def sql_bool(value: bool) -> str:
    return "TRUE" if bool(value) else "FALSE"


def connect_duckdb(cfg: dict, db_path: str = ":memory:") -> duckdb.DuckDBPyConnection:
    con = duckdb.connect(database=db_path)
    runtime = cfg.get("runtime", {})
    threads = int(runtime.get("n_jobs", 4))
    con.execute(f"PRAGMA threads={max(1, threads)}")
    memory_limit_mb = int(runtime.get("duckdb_memory_limit_mb", 4096))
    con.execute(f"PRAGMA memory_limit='{max(512, memory_limit_mb)}MB'")
    con.execute("PRAGMA preserve_insertion_order=FALSE")
    temp_dir = Path(runtime.get("duckdb_temp_dir", "data_interim/duckdb_tmp"))
    temp_dir.mkdir(parents=True, exist_ok=True)
    con.execute(f"PRAGMA temp_directory={sql_path(temp_dir)}")
    return con
