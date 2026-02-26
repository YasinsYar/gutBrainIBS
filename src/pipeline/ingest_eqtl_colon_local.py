from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd
import yaml

from src.utils.io import ensure_parent, sha256_file


def stage_file(src: Path, dst: Path) -> str:
    ensure_parent(dst)
    if dst.exists() and src.resolve() == dst.resolve():
        return "existing"
    if dst.exists():
        try:
            src_stat = src.stat()
            dst_stat = dst.stat()
            if dst_stat.st_size == src_stat.st_size and dst_stat.st_mtime_ns >= src_stat.st_mtime_ns:
                return "existing"
        except OSError:
            pass
    shutil.copy2(src, dst)
    return "copied"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--out-sigmoid", required=True)
    parser.add_argument("--out-transverse", required=True)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle)

    src_sig = Path(cfg["eqtl"]["colon_local"]["colon_sigmoid"])
    src_trv = Path(cfg["eqtl"]["colon_local"]["colon_transverse"])
    dst_sig = Path(args.out_sigmoid)
    dst_trv = Path(args.out_transverse)

    status_sig = stage_file(src_sig, dst_sig)
    status_trv = stage_file(src_trv, dst_trv)

    rows = [
        {
            "dataset_id": "QTD000226",
            "tissue": "colon_sigmoid",
            "source": "eqtl_catalogue",
            "local_path": str(dst_sig),
            "status": status_sig,
            "sha256": sha256_file(dst_sig),
        },
        {
            "dataset_id": "QTD000231",
            "tissue": "colon_transverse",
            "source": "eqtl_catalogue",
            "local_path": str(dst_trv),
            "status": status_trv,
            "sha256": sha256_file(dst_trv),
        },
    ]

    manifest = pd.DataFrame(rows)
    ensure_parent(args.manifest)
    manifest.to_csv(args.manifest, sep="\t", index=False)


if __name__ == "__main__":
    main()
