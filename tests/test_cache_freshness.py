from __future__ import annotations

import os
from pathlib import Path

from src.pipeline.ingest_eqtl_colon_local import stage_file
from src.pipeline.parse_eqtl_to_parquet import part_filename


def test_stage_file_recopies_when_source_is_newer_same_size(tmp_path: Path) -> None:
    src = tmp_path / "src.tsv.gz"
    dst = tmp_path / "dst.tsv.gz"
    src.write_bytes(b"AAAA")
    dst.write_bytes(b"BBBB")

    # keep size equal, make dst newer first -> no copy
    src_time = 1_700_000_000
    dst_time = 1_700_000_100
    src.touch()
    dst.touch()
    os.utime(src, (src_time, src_time))
    os.utime(dst, (dst_time, dst_time))
    status = stage_file(src, dst)
    assert status == "existing"
    assert dst.read_bytes() == b"BBBB"

    # now source becomes newer, stage_file must copy even with same size
    src.write_bytes(b"CCCC")
    newer = 1_700_000_300
    os.utime(src, (newer, newer))
    status = stage_file(src, dst)
    assert status == "copied"
    assert dst.read_bytes() == b"CCCC"


def test_part_filename_depends_on_source_fingerprint() -> None:
    base = {
        "dataset_id": "QTD000226",
        "tissue": "colon_sigmoid",
        "path": Path("data_raw/qtl/colon_sigmoid.all.tsv.gz"),
        "size_bytes": 1234,
    }
    a = part_filename({**base, "sha256": "abc"})
    b = part_filename({**base, "sha256": "def"})
    assert a != b
