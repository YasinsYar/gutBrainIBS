from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import h5py
import pandas as pd
import yaml

from src.utils.io import ensure_parent, sha256_file, write_json, write_tsv


def check_gzip(path: Path) -> tuple[bool, str]:
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as handle:
            _ = handle.readline()
        return True, "ok"
    except Exception as exc:
        return False, f"gzip_error: {exc}"


def check_h5ad(path: Path) -> tuple[bool, str]:
    try:
        with h5py.File(path, "r"):
            pass
        return True, "ok"
    except Exception as exc:
        return False, f"h5ad_error: {exc}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--inventory", required=True)
    parser.add_argument("--summary-json", required=True)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle)

    rows: list[dict] = []

    colon_inputs = cfg["eqtl"]["colon_local"]
    for tissue, raw_path in colon_inputs.items():
        path = Path(raw_path)
        exists = path.exists()
        ok, note = check_gzip(path) if exists else (False, "missing")
        rows.append(
            {
                "section": "eqtl_colon",
                "id": tissue,
                "path": str(path),
                "exists": exists,
                "readable": ok,
                "size_bytes": path.stat().st_size if exists else 0,
                "note": note,
            }
        )

    for idx, raw_path in enumerate(cfg["sc"]["inputs"]):
        path = Path(raw_path)
        exists = path.exists()
        ok, note = check_h5ad(path) if exists else (False, "missing")
        rows.append(
            {
                "section": "scrna",
                "id": f"sc_{idx + 1}",
                "path": str(path),
                "exists": exists,
                "readable": ok,
                "size_bytes": path.stat().st_size if exists else 0,
                "note": note,
            }
        )

    with open(cfg["gwas"]["manifest"], "r", encoding="utf-8") as handle:
        manifest = yaml.safe_load(handle)

    for entry in manifest["traits"]:
        path = Path(entry["file_path"])
        exists = path.exists()
        if exists:
            ok = True
            note = "ok"
        else:
            ok = not entry.get("is_mandatory", True)
            note = "missing_mandatory" if entry.get("is_mandatory", True) else "missing_optional"
        rows.append(
            {
                "section": "gwas",
                "id": entry["trait_id"],
                "path": str(path),
                "exists": exists,
                "readable": ok,
                "size_bytes": path.stat().st_size if exists else 0,
                "note": note,
            }
        )

    leads_path = Path("data_raw/ibsgwas_leads.tsv")
    rows.append(
        {
            "section": "gwas",
            "id": "ibs_leads",
            "path": str(leads_path),
            "exists": leads_path.exists(),
            "readable": leads_path.exists(),
            "size_bytes": leads_path.stat().st_size if leads_path.exists() else 0,
            "note": "ok" if leads_path.exists() else "missing",
        }
    )

    inv = pd.DataFrame(rows)
    ensure_parent(args.inventory)
    write_tsv(inv, args.inventory)

    n_missing_colon = int(((inv["section"] == "eqtl_colon") & (~inv["exists"])).sum())
    n_unreadable_colon = int(((inv["section"] == "eqtl_colon") & inv["exists"] & (~inv["readable"])).sum())
    n_missing_mandatory_gwas = int(((inv["section"] == "gwas") & (inv["note"] == "missing_mandatory")).sum())

    summary = {
        "n_items": int(len(inv)),
        "n_missing": int((~inv["exists"]).sum()),
        "n_unreadable": int((inv["exists"] & ~inv["readable"]).sum()),
        "n_missing_colon": n_missing_colon,
        "n_unreadable_colon": n_unreadable_colon,
        "n_missing_mandatory_gwas": n_missing_mandatory_gwas,
        "warnings": inv.loc[inv["note"] != "ok", ["section", "id", "note", "path"]].to_dict(orient="records"),
    }
    write_json(summary, args.summary_json)

    if n_missing_colon > 0 or n_unreadable_colon > 0 or n_missing_mandatory_gwas > 0:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
