from __future__ import annotations

import argparse
import copy
from pathlib import Path

import numpy as np
import pandas as pd
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


def p12_tag(value: float) -> str:
    return f"{value:.0e}".replace("+", "")


def add_baseline_overlap_stats(summary: pd.DataFrame, detail: pd.DataFrame, baseline_id: str) -> pd.DataFrame:
    if summary.empty:
        return summary
    out = summary.copy()
    out["n_overlap_with_baseline"] = 0
    out["pp4_spearman_vs_baseline"] = np.nan
    out["top_hit_changed_vs_baseline"] = False

    base = detail.loc[detail["scenario_id"] == baseline_id, ["locus_id", "tissue", "gene_id", "PP4"]].copy()
    if base.empty:
        return out
    base = base.rename(columns={"PP4": "PP4_baseline"})
    base_top = base.sort_values("PP4_baseline", ascending=False).head(1)
    if base_top.empty:
        base_key = ("NA", "NA", "NA")
    else:
        r = base_top.iloc[0]
        base_key = (str(r["locus_id"]), str(r["tissue"]), str(r["gene_id"]))

    for idx, row in out.iterrows():
        scenario_id = str(row["scenario_id"])
        sub = detail.loc[detail["scenario_id"] == scenario_id, ["locus_id", "tissue", "gene_id", "PP4"]].copy()
        if sub.empty:
            continue
        merged = sub.merge(base, on=["locus_id", "tissue", "gene_id"], how="inner")
        out.at[idx, "n_overlap_with_baseline"] = int(len(merged))
        if len(merged) >= 3:
            out.at[idx, "pp4_spearman_vs_baseline"] = float(merged["PP4"].corr(merged["PP4_baseline"], method="spearman"))
        top = sub.sort_values("PP4", ascending=False).head(1)
        if not top.empty:
            t = top.iloc[0]
            key = (str(t["locus_id"]), str(t["tissue"]), str(t["gene_id"]))
            out.at[idx, "top_hit_changed_vs_baseline"] = key != base_key
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--eqtl", required=True)
    parser.add_argument("--gwas", required=True)
    parser.add_argument("--out-summary", required=True)
    parser.add_argument("--out-brain-detail", required=True)
    args = parser.parse_args()

    cfg = yaml.safe_load(open(args.config, "r", encoding="utf-8"))
    sensitivity_cfg = cfg.get("sensitivity", {})
    coloc_cfg = sensitivity_cfg.get("coloc", {})

    windows_bp = sorted({int(x) for x in coloc_cfg.get("windows_bp", [250000, 500000, 1000000])})
    p12_values = sorted({float(x) for x in coloc_cfg.get("p12_values", [1.0e-6, 1.0e-5, 1.0e-4])})
    mhc_start = int(coloc_cfg.get("mhc_start", 25000000))
    mhc_end = int(coloc_cfg.get("mhc_end", 34000000))
    baseline_window = int(cfg["analysis"]["window_bp"])
    baseline_p12 = float(cfg["analysis"]["coloc_priors"]["p12"])
    progress_every = int(cfg.get("runtime", {}).get("progress_every_row_groups", 20))

    ibs_trait_id = str(cfg["gwas"]["ibs_trait_id"])
    leads = load_leads(Path("data_raw/ibsgwas_leads.tsv"), Path(args.gwas), ibs_trait_id, progress_every)
    brain_tissues = scan_brain_tissues(Path(args.eqtl))

    if leads.empty or not brain_tissues:
        pd.DataFrame(
            [
                {
                    "scenario_id": "baseline_all",
                    "window_bp": baseline_window,
                    "p12": baseline_p12,
                    "mhc_mode": "all",
                    "n_rows": 0,
                    "n_pp4_ge_0_7": 0,
                    "n_pp4_ge_0_5": 0,
                    "pp4_mean": np.nan,
                    "pp4_median": np.nan,
                    "pp4_max": np.nan,
                    "top_locus_id": "NA",
                    "top_gene_id": "NA",
                    "top_tissue": "NA",
                    "status": "no_leads_or_no_brain_tissues",
                }
            ]
        ).to_csv(args.out_summary, sep="\t", index=False)
        pd.DataFrame(columns=["scenario_id", "locus_id", "lead_snp", "tissue", "gene_id", "PP4"]).to_csv(
            args.out_brain_detail, sep="\t", index=False
        )
        return

    max_window = max(windows_bp + [baseline_window])
    windows = build_windows(leads, max_window)
    eqtl_local = collect_eqtl_in_windows(Path(args.eqtl), set(brain_tissues), windows, progress_every)
    if eqtl_local.empty:
        pd.DataFrame(
            [
                {
                    "scenario_id": "baseline_all",
                    "window_bp": baseline_window,
                    "p12": baseline_p12,
                    "mhc_mode": "all",
                    "n_rows": 0,
                    "n_pp4_ge_0_7": 0,
                    "n_pp4_ge_0_5": 0,
                    "pp4_mean": np.nan,
                    "pp4_median": np.nan,
                    "pp4_max": np.nan,
                    "top_locus_id": "NA",
                    "top_gene_id": "NA",
                    "top_tissue": "NA",
                    "status": "no_eqtl_in_windows",
                }
            ]
        ).to_csv(args.out_summary, sep="\t", index=False)
        pd.DataFrame(columns=["scenario_id", "locus_id", "lead_snp", "tissue", "gene_id", "PP4"]).to_csv(
            args.out_brain_detail, sep="\t", index=False
        )
        return

    lead_pos = eqtl_local["locus_id"].astype(str).str.split(":", n=1).str[1]
    eqtl_local["lead_pos"] = pd.to_numeric(lead_pos, errors="coerce").astype("Int64")
    eqtl_local = eqtl_local.dropna(subset=["lead_pos"]).copy()
    eqtl_local["lead_pos"] = eqtl_local["lead_pos"].astype(int)
    eqtl_local["distance_bp"] = (eqtl_local["pos"].astype(int) - eqtl_local["lead_pos"]).abs()

    keys_by_chr = build_keys_by_chr(eqtl_local)
    gwas_map = collect_gwas_map(Path(args.gwas), ibs_trait_id, keys_by_chr, progress_every)

    scenarios: list[dict] = []
    for window_bp in windows_bp:
        for p12 in p12_values:
            scenarios.append(
                {
                    "scenario_id": f"w{window_bp}_p12_{p12_tag(p12)}_all",
                    "window_bp": int(window_bp),
                    "p12": float(p12),
                    "mhc_mode": "all",
                }
            )
    scenarios.extend(
        [
            {
                "scenario_id": f"w{baseline_window}_p12_{p12_tag(baseline_p12)}_exclude_mhc",
                "window_bp": baseline_window,
                "p12": baseline_p12,
                "mhc_mode": "exclude_mhc",
            },
            {
                "scenario_id": f"w{baseline_window}_p12_{p12_tag(baseline_p12)}_chr6_only",
                "window_bp": baseline_window,
                "p12": baseline_p12,
                "mhc_mode": "chr6_only",
            },
        ]
    )
    unique: dict[str, dict] = {}
    for sc in scenarios:
        unique[sc["scenario_id"]] = sc
    scenarios = list(unique.values())

    summary_rows: list[dict] = []
    detail_frames: list[pd.DataFrame] = []

    for idx, scenario in enumerate(scenarios, start=1):
        sc_id = str(scenario["scenario_id"])
        window_bp = int(scenario["window_bp"])
        p12 = float(scenario["p12"])
        mhc_mode = str(scenario["mhc_mode"])

        local = eqtl_local.loc[eqtl_local["distance_bp"] <= window_bp].copy()
        if mhc_mode == "exclude_mhc":
            local = local.loc[~((local["chr"] == "chr6") & local["pos"].between(mhc_start, mhc_end))]
        elif mhc_mode == "chr6_only":
            local = local.loc[local["chr"] == "chr6"]

        if local.empty:
            summary_rows.append(
                {
                    "scenario_id": sc_id,
                    "window_bp": window_bp,
                    "p12": p12,
                    "mhc_mode": mhc_mode,
                    "n_rows": 0,
                    "n_pp4_ge_0_7": 0,
                    "n_pp4_ge_0_5": 0,
                    "pp4_mean": np.nan,
                    "pp4_median": np.nan,
                    "pp4_max": np.nan,
                    "top_locus_id": "NA",
                    "top_gene_id": "NA",
                    "top_tissue": "NA",
                    "status": "no_rows_after_filter",
                }
            )
            continue

        local_cfg = copy.deepcopy(cfg)
        local_cfg["analysis"]["coloc_priors"]["p12"] = p12
        table = build_table(local, gwas_map, brain_tissues, local_cfg)

        if table.empty:
            summary_rows.append(
                {
                    "scenario_id": sc_id,
                    "window_bp": window_bp,
                    "p12": p12,
                    "mhc_mode": mhc_mode,
                    "n_rows": 0,
                    "n_pp4_ge_0_7": 0,
                    "n_pp4_ge_0_5": 0,
                    "pp4_mean": np.nan,
                    "pp4_median": np.nan,
                    "pp4_max": np.nan,
                    "top_locus_id": "NA",
                    "top_gene_id": "NA",
                    "top_tissue": "NA",
                    "status": "no_coloc_rows",
                }
            )
            continue

        table["scenario_id"] = sc_id
        table["window_bp"] = window_bp
        table["p12"] = p12
        table["mhc_mode"] = mhc_mode
        detail_frames.append(table)

        top = table.sort_values("PP4", ascending=False).iloc[0]
        summary_rows.append(
            {
                "scenario_id": sc_id,
                "window_bp": window_bp,
                "p12": p12,
                "mhc_mode": mhc_mode,
                "n_rows": int(len(table)),
                "n_pp4_ge_0_7": int((table["PP4"] >= 0.7).sum()),
                "n_pp4_ge_0_5": int((table["PP4"] >= 0.5).sum()),
                "pp4_mean": float(table["PP4"].mean()),
                "pp4_median": float(table["PP4"].median()),
                "pp4_max": float(table["PP4"].max()),
                "top_locus_id": str(top["locus_id"]),
                "top_gene_id": str(top["gene_id"]),
                "top_tissue": str(top["tissue"]),
                "status": "ok",
            }
        )
        print(f"[coloc_sensitivity] scenario={idx}/{len(scenarios)} id={sc_id} rows={len(table)}", flush=True)

    summary = pd.DataFrame(summary_rows)
    if detail_frames:
        detail = pd.concat(detail_frames, ignore_index=True)
    else:
        detail = pd.DataFrame(columns=["scenario_id", "locus_id", "lead_snp", "tissue", "gene_id", "PP4"])

    baseline_id = f"w{baseline_window}_p12_{p12_tag(baseline_p12)}_all"
    summary = add_baseline_overlap_stats(summary, detail, baseline_id)
    summary = summary.sort_values(["mhc_mode", "window_bp", "p12"]).reset_index(drop=True)

    Path(args.out_summary).parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.out_summary, sep="\t", index=False)
    detail.to_csv(args.out_brain_detail, sep="\t", index=False)


if __name__ == "__main__":
    main()
