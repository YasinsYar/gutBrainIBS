from __future__ import annotations

import argparse
import copy
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import yaml

from src.pipeline.coloc_colon_and_brain import (
    build_keys_by_chr,
    build_table,
    build_windows,
    collect_eqtl_in_windows,
    collect_gwas_map,
    load_leads,
    scan_brain_tissues,
)
from src.utils.genetics import is_palindromic


def label_from_criteria(crit_rg: bool, crit_coloc: bool, crit_mr: bool) -> str:
    if crit_rg and crit_coloc and crit_mr:
        return "Supported gut-brain genetic contribution"
    if crit_rg or crit_coloc or crit_mr:
        return "Partial support"
    return "Insufficient support"


def collect_palindromic_sig_in_windows(
    eqtl_parsed: Path,
    tissues: set[str],
    windows: dict[str, list[dict]],
    p_thr: float,
    progress_every: int,
) -> pd.DataFrame:
    cols = ["chr", "pos", "variant", "gene_id", "tissue", "beta", "se", "pvalue", "ref", "alt"]
    pf = pq.ParquetFile(eqtl_parsed)
    out_chunks: list[pd.DataFrame] = []
    valid_chr = set(windows.keys())
    started = time.time()
    for rg in range(pf.metadata.num_row_groups):
        chunk = pf.read_row_group(rg, columns=cols).to_pandas()
        chunk = chunk.loc[chunk["tissue"].isin(tissues) & chunk["chr"].isin(valid_chr)].copy()
        if chunk.empty:
            if (rg + 1) % progress_every == 0:
                elapsed = time.time() - started
                print(
                    f"[harmonization_sensitivity] row_groups={rg + 1}/{pf.metadata.num_row_groups} matched_chunks={len(out_chunks)} elapsed={elapsed:.1f}s",
                    flush=True,
                )
            continue
        chunk["pos"] = pd.to_numeric(chunk["pos"], errors="coerce")
        chunk["beta"] = pd.to_numeric(chunk["beta"], errors="coerce")
        chunk["se"] = pd.to_numeric(chunk["se"], errors="coerce")
        chunk["pvalue"] = pd.to_numeric(chunk["pvalue"], errors="coerce")
        chunk = chunk.dropna(subset=["chr", "pos", "variant", "gene_id", "tissue", "beta", "se", "pvalue", "ref", "alt"])
        if chunk.empty:
            if (rg + 1) % progress_every == 0:
                elapsed = time.time() - started
                print(
                    f"[harmonization_sensitivity] row_groups={rg + 1}/{pf.metadata.num_row_groups} matched_chunks={len(out_chunks)} elapsed={elapsed:.1f}s",
                    flush=True,
                )
            continue
        chunk = chunk.loc[chunk["pvalue"] < p_thr].copy()
        if chunk.empty:
            if (rg + 1) % progress_every == 0:
                elapsed = time.time() - started
                print(
                    f"[harmonization_sensitivity] row_groups={rg + 1}/{pf.metadata.num_row_groups} matched_chunks={len(out_chunks)} elapsed={elapsed:.1f}s",
                    flush=True,
                )
            continue
        chunk["pos"] = chunk["pos"].astype(int)
        chunk["ref"] = chunk["ref"].astype(str).str.upper()
        chunk["alt"] = chunk["alt"].astype(str).str.upper()
        pal_mask = [is_palindromic(r, a) for r, a in zip(chunk["ref"], chunk["alt"])]
        chunk = chunk.loc[pal_mask].copy()
        if chunk.empty:
            if (rg + 1) % progress_every == 0:
                elapsed = time.time() - started
                print(
                    f"[harmonization_sensitivity] row_groups={rg + 1}/{pf.metadata.num_row_groups} matched_chunks={len(out_chunks)} elapsed={elapsed:.1f}s",
                    flush=True,
                )
            continue
        for chrom, sub in chunk.groupby("chr"):
            for window in windows.get(str(chrom), []):
                local = sub.loc[sub["pos"].between(window["start"], window["end"])].copy()
                if local.empty:
                    continue
                local["locus_id"] = window["locus_id"]
                local["lead_snp"] = window["lead_snp"]
                local["lead_beta"] = window["lead_beta"]
                local["lead_se"] = window["lead_se"]
                out_chunks.append(local)
        if (rg + 1) % progress_every == 0:
            elapsed = time.time() - started
            print(
                f"[harmonization_sensitivity] row_groups={rg + 1}/{pf.metadata.num_row_groups} matched_chunks={len(out_chunks)} elapsed={elapsed:.1f}s",
                flush=True,
            )
    if not out_chunks:
        return pd.DataFrame(
            columns=[
                "chr",
                "pos",
                "variant",
                "gene_id",
                "tissue",
                "beta",
                "se",
                "pvalue",
                "locus_id",
                "lead_snp",
                "lead_beta",
                "lead_se",
            ]
        )
    out = pd.concat(out_chunks, ignore_index=True)
    out = out.sort_values("pvalue").drop_duplicates(
        subset=["locus_id", "lead_snp", "tissue", "gene_id", "variant", "chr", "pos"],
        keep="first",
    )
    return out


def top_hit(frame: pd.DataFrame) -> tuple[str, str, str, float]:
    if frame.empty:
        return "NA", "NA", "NA", np.nan
    sub = frame.sort_values("PP4", ascending=False).iloc[0]
    return str(sub.get("locus_id", "NA")), str(sub.get("gene_id", "NA")), str(sub.get("tissue", "NA")), float(sub.get("PP4", np.nan))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--eqtl-parsed", required=True)
    parser.add_argument("--eqtl-significant", required=True)
    parser.add_argument("--gwas", required=True)
    parser.add_argument("--strict-coloc-brain", required=True)
    parser.add_argument("--strict-scorecard", default="")
    parser.add_argument("--out-table", required=True)
    parser.add_argument("--out-coloc-brain-nonstrict", required=True)
    args = parser.parse_args()

    cfg = yaml.safe_load(open(args.config, "r", encoding="utf-8"))
    ibs_trait_id = str(cfg["gwas"]["ibs_trait_id"])
    p_thr = float(cfg["analysis"]["p_threshold_eqtl"])
    window_bp = int(cfg["analysis"]["window_bp"])
    progress_every = int(cfg.get("runtime", {}).get("progress_every_row_groups", 20))

    leads = load_leads(Path("data_raw/ibsgwas_leads.tsv"), Path(args.gwas), ibs_trait_id, progress_every)
    brain_tissues = scan_brain_tissues(Path(args.eqtl_significant))
    if leads.empty or not brain_tissues:
        pd.DataFrame(
            [
                {
                    "status": "no_leads_or_brain_tissues",
                    "strict_n_rows": 0,
                    "nonstrict_n_rows": 0,
                    "n_added_palindromic_window_rows": 0,
                    "mr_instruments_changed": False,
                    "mr_note": "mr_not_affected_by_eqtl_palindromic_mode",
                }
            ]
        ).to_csv(args.out_table, sep="\t", index=False)
        pd.DataFrame(columns=["locus_id", "lead_snp", "tissue", "gene_id", "PP4"]).to_csv(
            args.out_coloc_brain_nonstrict, sep="\t", index=False
        )
        return

    windows = build_windows(leads, window_bp)
    strict_local = collect_eqtl_in_windows(Path(args.eqtl_significant), set(brain_tissues), windows, progress_every)
    extra_pal = collect_palindromic_sig_in_windows(Path(args.eqtl_parsed), set(brain_tissues), windows, p_thr, progress_every)
    nonstrict_local = strict_local.copy()
    if not extra_pal.empty:
        nonstrict_local = pd.concat([nonstrict_local, extra_pal], ignore_index=True)
        nonstrict_local = nonstrict_local.sort_values("pvalue").drop_duplicates(
            subset=["locus_id", "lead_snp", "tissue", "gene_id", "variant", "chr", "pos"],
            keep="first",
        )

    keys_by_chr = build_keys_by_chr(nonstrict_local if not nonstrict_local.empty else strict_local)
    gwas_map = collect_gwas_map(Path(args.gwas), ibs_trait_id, keys_by_chr, progress_every)

    if Path(args.strict_coloc_brain).exists():
        strict_coloc = pd.read_csv(args.strict_coloc_brain, sep="\t")
    else:
        strict_coloc = pd.DataFrame()
    if strict_coloc.empty:
        strict_coloc = build_table(strict_local, gwas_map, brain_tissues, cfg)
    strict_coloc = strict_coloc.copy()
    nonstrict_cfg = copy.deepcopy(cfg)
    nonstrict_cfg["analysis"]["exclude_palindromic"] = False
    nonstrict_coloc = build_table(nonstrict_local, gwas_map, brain_tissues, nonstrict_cfg)

    strict_top_locus, strict_top_gene, strict_top_tissue, strict_top_pp4 = top_hit(strict_coloc)
    nonstrict_top_locus, nonstrict_top_gene, nonstrict_top_tissue, nonstrict_top_pp4 = top_hit(nonstrict_coloc)
    top_changed = (strict_top_locus, strict_top_gene, strict_top_tissue) != (
        nonstrict_top_locus,
        nonstrict_top_gene,
        nonstrict_top_tissue,
    )

    crit_rg = False
    crit_mr = False
    strict_label = "NA"
    nonstrict_label = "NA"
    if args.strict_scorecard and Path(args.strict_scorecard).exists():
        scorecard = json.loads(Path(args.strict_scorecard).read_text(encoding="utf-8"))
        crit_rg = bool(scorecard.get("criteria", {}).get("genetic_correlation_supported", False))
        crit_mr = bool(scorecard.get("criteria", {}).get("mr_directional_support", False))
        strict_label = str(scorecard.get("interpretation_label", "NA"))

    strict_coloc_ok = bool((strict_coloc.get("PP4", pd.Series(dtype=float)) >= 0.7).any()) if not strict_coloc.empty else False
    nonstrict_coloc_ok = bool((nonstrict_coloc.get("PP4", pd.Series(dtype=float)) >= 0.7).any()) if not nonstrict_coloc.empty else False
    if strict_label == "NA":
        strict_label = label_from_criteria(crit_rg, strict_coloc_ok, crit_mr)
    nonstrict_label = label_from_criteria(crit_rg, nonstrict_coloc_ok, crit_mr)

    out_row = {
        "status": "ok",
        "strict_n_rows": int(len(strict_coloc)),
        "nonstrict_n_rows": int(len(nonstrict_coloc)),
        "strict_n_pp4_ge_0_7": int((strict_coloc["PP4"] >= 0.7).sum()) if not strict_coloc.empty else 0,
        "nonstrict_n_pp4_ge_0_7": int((nonstrict_coloc["PP4"] >= 0.7).sum()) if not nonstrict_coloc.empty else 0,
        "strict_pp4_max": float(strict_coloc["PP4"].max()) if not strict_coloc.empty else np.nan,
        "nonstrict_pp4_max": float(nonstrict_coloc["PP4"].max()) if not nonstrict_coloc.empty else np.nan,
        "strict_top_locus_id": strict_top_locus,
        "strict_top_gene_id": strict_top_gene,
        "strict_top_tissue": strict_top_tissue,
        "strict_top_pp4": strict_top_pp4,
        "nonstrict_top_locus_id": nonstrict_top_locus,
        "nonstrict_top_gene_id": nonstrict_top_gene,
        "nonstrict_top_tissue": nonstrict_top_tissue,
        "nonstrict_top_pp4": nonstrict_top_pp4,
        "top_hit_changed": bool(top_changed),
        "n_added_palindromic_window_rows": int(max(len(nonstrict_local) - len(strict_local), 0)),
        "strict_label": strict_label,
        "nonstrict_label": nonstrict_label,
        "mr_instruments_changed": False,
        "mr_note": "mr_not_affected_by_eqtl_palindromic_mode",
    }

    Path(args.out_table).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([out_row]).to_csv(args.out_table, sep="\t", index=False)
    nonstrict_coloc.to_csv(args.out_coloc_brain_nonstrict, sep="\t", index=False)


if __name__ == "__main__":
    main()
