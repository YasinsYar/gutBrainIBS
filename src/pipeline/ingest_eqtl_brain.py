from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd
import requests
import yaml

from src.utils.io import ensure_parent, sha256_file


def ftp_to_https(url: str) -> str:
    if url.startswith("ftp://ftp.ebi.ac.uk"):
        return url.replace("ftp://", "https://", 1)
    return url


def remote_size(url: str) -> int | None:
    try:
        with requests.get(url, stream=True, timeout=120) as response:
            response.raise_for_status()
            value = response.headers.get("content-length")
            if value is None:
                return None
            return int(value)
    except Exception:
        return None


def download_file(url: str, out_path: Path) -> str:
    ensure_parent(out_path)
    expected_size = remote_size(url)
    if out_path.exists() and out_path.stat().st_size > 0:
        local_size = out_path.stat().st_size
        if expected_size is None or local_size == expected_size:
            return "existing"
        out_path.unlink()
    with requests.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        with open(out_path, "wb") as handle:
            shutil.copyfileobj(response.raw, handle)
    if expected_size is not None and out_path.stat().st_size != expected_size:
        raise RuntimeError(
            f"size_mismatch_after_download local={out_path.stat().st_size} expected={expected_size}"
        )
    return "downloaded"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle)

    out_rows: list[dict] = []
    brain_cfg = cfg["eqtl"]["brain"]
    if not brain_cfg.get("enabled", True):
        pd.DataFrame(columns=["dataset_id", "tissue", "source", "local_path", "status", "sha256"]).to_csv(
            args.manifest, sep="\t", index=False
        )
        return

    metadata = pd.read_csv(brain_cfg["source_table_url"], sep="\t")
    subset = metadata[
        (metadata["study_label"].str.lower() == "gtex")
        & (metadata["quant_method"].str.lower() == brain_cfg.get("quant_method", "ge"))
    ].copy()

    keywords = [x.lower() for x in brain_cfg.get("tissue_keywords", ["brain"])]
    mask = subset["tissue_label"].str.lower().apply(lambda x: any(k in x for k in keywords))
    subset = subset[mask].sort_values("dataset_id").head(int(brain_cfg.get("max_datasets", 4)))
    if subset.empty:
        raise RuntimeError("no_brain_datasets_selected_from_metadata")

    out_dir = Path("data_raw/qtl/brain")
    out_dir.mkdir(parents=True, exist_ok=True)

    for row in subset.itertuples(index=False):
        dataset_id = row.dataset_id
        tissue = str(row.sample_group)
        source_url = ftp_to_https(str(row.ftp_path))
        local_path = out_dir / f"{dataset_id}.all.tsv.gz"
        status = "failed"
        sha = ""
        try:
            status = download_file(source_url, local_path)
            sha = sha256_file(local_path)
        except Exception as exc:
            status = f"failed: {exc}"
        out_rows.append(
            {
                "dataset_id": dataset_id,
                "tissue": tissue,
                "source": "eqtl_catalogue",
                "local_path": str(local_path),
                "status": status,
                "sha256": sha,
                "url": source_url,
            }
        )

    manifest = pd.DataFrame(out_rows)
    manifest.to_csv(args.manifest, sep="\t", index=False)
    n_ok = int(manifest["status"].astype(str).isin(["existing", "downloaded"]).sum()) if not manifest.empty else 0
    if n_ok == 0:
        raise RuntimeError("brain_eqtl_ingest_failed_all_datasets")


if __name__ == "__main__":
    main()
