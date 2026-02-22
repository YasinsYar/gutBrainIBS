from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import yaml

from src.utils.genetics import clump_by_distance
from src.utils.stats import egger_mr, ivw_mr, weighted_median


def load_trait_instruments(parquet_path: str, trait_id: str, p_thr: float, progress_every: int) -> pd.DataFrame:
    pf = pq.ParquetFile(parquet_path)
    best: dict[tuple[str, int], tuple[float, float, float]] = {}
    cols = ["trait_id", "chr", "pos", "beta", "se", "p"]
    for rg in range(pf.metadata.num_row_groups):
        chunk = pf.read_row_group(rg, columns=cols).to_pandas()
        chunk = chunk.loc[(chunk["trait_id"] == trait_id)].copy()
        if chunk.empty:
            continue
        chunk["pos"] = pd.to_numeric(chunk["pos"], errors="coerce")
        chunk["beta"] = pd.to_numeric(chunk["beta"], errors="coerce")
        chunk["se"] = pd.to_numeric(chunk["se"], errors="coerce")
        chunk["p"] = pd.to_numeric(chunk["p"], errors="coerce")
        chunk = chunk.dropna(subset=["chr", "pos", "beta", "se", "p"])
        if chunk.empty:
            continue
        chunk = chunk.loc[chunk["p"] < p_thr]
        if chunk.empty:
            continue
        chunk["chr"] = chunk["chr"].astype(str)
        chunk["pos"] = chunk["pos"].astype(int)
        chunk = chunk.sort_values("p").drop_duplicates(subset=["chr", "pos"], keep="first")
        for row in chunk.itertuples(index=False):
            key = (str(row.chr), int(row.pos))
            pvalue = float(row.p)
            prev = best.get(key)
            if prev is None or pvalue < prev[2]:
                best[key] = (float(row.beta), float(row.se), pvalue)
        if (rg + 1) % progress_every == 0:
            print(
                f"[bidirectional_mr] instruments trait={trait_id} row_groups={rg + 1}/{pf.metadata.num_row_groups} selected={len(best)}",
                flush=True,
            )
    if not best:
        return pd.DataFrame(columns=["chr", "pos", "beta", "se", "p"])
    rows = [{"chr": k[0], "pos": k[1], "beta": v[0], "se": v[1], "p": v[2]} for k, v in best.items()]
    out = pd.DataFrame(rows).sort_values("p")
    return out


def positions_by_chr(frame: pd.DataFrame) -> dict[str, set[int]]:
    if frame.empty:
        return {}
    out: dict[str, set[int]] = {}
    for chrom, sub in frame.groupby("chr"):
        out[str(chrom)] = set(sub["pos"].astype(int).tolist())
    return out


def load_outcome_for_positions(parquet_path: str, trait_id: str, pos_by_chr: dict[str, set[int]], progress_every: int) -> pd.DataFrame:
    if not pos_by_chr:
        return pd.DataFrame(columns=["chr", "pos", "beta", "se", "p"])
    pf = pq.ParquetFile(parquet_path)
    best: dict[tuple[str, int], tuple[float, float, float]] = {}
    cols = ["trait_id", "chr", "pos", "beta", "se", "p"]
    valid_chr = set(pos_by_chr.keys())
    for rg in range(pf.metadata.num_row_groups):
        chunk = pf.read_row_group(rg, columns=cols).to_pandas()
        chunk = chunk.loc[(chunk["trait_id"] == trait_id) & (chunk["chr"].isin(valid_chr))].copy()
        if chunk.empty:
            continue
        chunk["pos"] = pd.to_numeric(chunk["pos"], errors="coerce")
        chunk["beta"] = pd.to_numeric(chunk["beta"], errors="coerce")
        chunk["se"] = pd.to_numeric(chunk["se"], errors="coerce")
        chunk["p"] = pd.to_numeric(chunk["p"], errors="coerce")
        chunk = chunk.dropna(subset=["chr", "pos", "beta", "se", "p"])
        if chunk.empty:
            continue
        chunk["chr"] = chunk["chr"].astype(str)
        chunk["pos"] = chunk["pos"].astype(int)
        for chrom, sub in chunk.groupby("chr"):
            allowed = pos_by_chr.get(str(chrom), set())
            sub = sub.loc[sub["pos"].isin(allowed)]
            if sub.empty:
                continue
            sub = sub.sort_values("p").drop_duplicates(subset=["chr", "pos"], keep="first")
            for row in sub.itertuples(index=False):
                key = (str(row.chr), int(row.pos))
                pvalue = float(row.p)
                prev = best.get(key)
                if prev is None or pvalue < prev[2]:
                    best[key] = (float(row.beta), float(row.se), pvalue)
        if (rg + 1) % progress_every == 0:
            print(
                f"[bidirectional_mr] outcome trait={trait_id} row_groups={rg + 1}/{pf.metadata.num_row_groups} matched={len(best)}",
                flush=True,
            )
    if not best:
        return pd.DataFrame(columns=["chr", "pos", "beta", "se", "p"])
    rows = [{"chr": k[0], "pos": k[1], "beta": v[0], "se": v[1], "p": v[2]} for k, v in best.items()]
    return pd.DataFrame(rows)


def run_direction(exposure_sig: pd.DataFrame, outcome: pd.DataFrame, direction: str, kb: int, alpha: float) -> dict:
    if exposure_sig.empty:
        return {
            "direction": direction,
            "n_instruments": 0,
            "ivw_beta": np.nan,
            "ivw_se": np.nan,
            "ivw_p": np.nan,
            "wm_beta": np.nan,
            "egger_beta": np.nan,
            "egger_p": np.nan,
            "egger_intercept": np.nan,
            "egger_intercept_p": np.nan,
            "consistent_direction": False,
            "pleiotropy_flag": False,
            "status": "no_instruments",
        }

    exp_clumped = clump_by_distance(exposure_sig, p_col="p", chr_col="chr", pos_col="pos", kb=kb)
    merged = exp_clumped.merge(outcome, on=["chr", "pos"], suffixes=("_exp", "_out"))
    merged = merged.dropna(subset=["beta_exp", "se_exp", "beta_out", "se_out"])

    if len(merged) < 3:
        return {
            "direction": direction,
            "n_instruments": int(len(merged)),
            "ivw_beta": np.nan,
            "ivw_se": np.nan,
            "ivw_p": np.nan,
            "wm_beta": np.nan,
            "egger_beta": np.nan,
            "egger_p": np.nan,
            "egger_intercept": np.nan,
            "egger_intercept_p": np.nan,
            "consistent_direction": False,
            "pleiotropy_flag": False,
            "status": "insufficient_overlap",
        }

    beta_exp = merged["beta_exp"].to_numpy()
    se_exp = merged["se_exp"].to_numpy()
    beta_out = merged["beta_out"].to_numpy()
    se_out = merged["se_out"].to_numpy()

    ivw_beta, ivw_se, ivw_p = ivw_mr(beta_exp, se_exp, beta_out, se_out)
    ratio = beta_out / beta_exp
    var_ratio = np.square(se_out / beta_exp) + np.square(beta_out * se_exp / np.square(beta_exp))
    w = 1.0 / np.clip(var_ratio, 1e-12, None)
    wm_beta = weighted_median(ratio, w)

    egger_beta, egger_p, egger_intercept, egger_intercept_p = egger_mr(beta_exp, se_exp, beta_out, se_out)
    consistent = (np.sign(ivw_beta) == np.sign(wm_beta)) and (ivw_p < alpha)

    return {
        "direction": direction,
        "n_instruments": int(len(merged)),
        "ivw_beta": float(ivw_beta),
        "ivw_se": float(ivw_se),
        "ivw_p": float(ivw_p),
        "wm_beta": float(wm_beta),
        "egger_beta": float(egger_beta),
        "egger_p": float(egger_p),
        "egger_intercept": float(egger_intercept),
        "egger_intercept_p": float(egger_intercept_p),
        "consistent_direction": bool(consistent),
        "pleiotropy_flag": bool(egger_intercept_p < 0.05),
        "status": "ok",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--in-parquet", required=True)
    parser.add_argument("--out-table", required=True)
    args = parser.parse_args()

    cfg = yaml.safe_load(open(args.config, "r", encoding="utf-8"))
    traits_cfg = yaml.safe_load(open(cfg["gwas"]["traits"], "r", encoding="utf-8"))["traits"]

    ibs_id = cfg["gwas"]["ibs_trait_id"]
    p_thr = float(cfg["analysis"]["mr_instrument_p"])
    kb = int(cfg["analysis"]["mr_clump_kb"])
    alpha = float(cfg["analysis"]["fdr_alpha"])
    progress_every = int(cfg.get("runtime", {}).get("progress_every_row_groups", 20))

    rows = []
    psych_traits = [t for t in traits_cfg if t.get("group") == "psychiatric"]
    ibs_sig = load_trait_instruments(args.in_parquet, ibs_id, p_thr, progress_every)

    for trait in psych_traits:
        t_id = trait["trait_id"]
        psych_sig = load_trait_instruments(args.in_parquet, t_id, p_thr, progress_every)

        ibs_for_psych = load_outcome_for_positions(args.in_parquet, ibs_id, positions_by_chr(psych_sig), progress_every)
        psych_for_ibs = load_outcome_for_positions(args.in_parquet, t_id, positions_by_chr(ibs_sig), progress_every)

        r1 = run_direction(psych_sig, ibs_for_psych, f"{t_id}->IBS", kb, alpha)
        r1.update({"exposure_trait": t_id, "outcome_trait": ibs_id})
        rows.append(r1)

        r2 = run_direction(ibs_sig, psych_for_ibs, f"IBS->{t_id}", kb, alpha)
        r2.update({"exposure_trait": ibs_id, "outcome_trait": t_id})
        rows.append(r2)

    out = pd.DataFrame(rows)
    out.to_csv(args.out_table, sep="\t", index=False)


if __name__ == "__main__":
    main()
