from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
import yaml
from scipy.stats import chi2

from src.pipeline.bidirectional_mr import load_outcome_for_positions, load_trait_instruments, positions_by_chr
from src.utils.genetics import clump_by_distance
from src.utils.stats import ivw_mr


def prepare_merged(exposure_sig: pd.DataFrame, outcome: pd.DataFrame, kb: int) -> pd.DataFrame:
    if exposure_sig.empty or outcome.empty:
        return pd.DataFrame()
    exp_clumped = clump_by_distance(exposure_sig, p_col="p", chr_col="chr", pos_col="pos", kb=kb)
    merged = exp_clumped.merge(outcome, on=["chr", "pos"], suffixes=("_exp", "_out"))
    merged = merged.dropna(subset=["beta_exp", "se_exp", "beta_out", "se_out"]).copy()
    if merged.empty:
        return merged
    merged = merged.loc[merged["beta_exp"].abs() > 1e-12].copy()
    return merged


def cochran_q(beta_exp: np.ndarray, se_exp: np.ndarray, beta_out: np.ndarray, se_out: np.ndarray) -> tuple[float, int, float]:
    ratio = beta_out / beta_exp
    var_ratio = np.square(se_out / beta_exp) + np.square(beta_out * se_exp / np.square(beta_exp))
    w = 1.0 / np.clip(var_ratio, 1e-12, None)
    theta_ivw = float(np.sum(w * ratio) / np.sum(w))
    q = float(np.sum(w * np.square(ratio - theta_ivw)))
    dof = max(len(ratio) - 1, 0)
    p = float(chi2.sf(q, dof)) if dof > 0 else np.nan
    return q, dof, p


def leave_one_out_ivw(beta_exp: np.ndarray, se_exp: np.ndarray, beta_out: np.ndarray, se_out: np.ndarray, full_ivw: float) -> tuple[float, int]:
    n = len(beta_exp)
    if n < 4:
        return np.nan, 0
    betas = []
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        b, _, _ = ivw_mr(beta_exp[mask], se_exp[mask], beta_out[mask], se_out[mask])
        betas.append(float(b))
    loo = np.asarray(betas, dtype=float)
    max_delta = float(np.max(np.abs(loo - full_ivw))) if len(loo) else np.nan
    sign_flip = int(np.sum(np.sign(loo) != np.sign(full_ivw))) if np.isfinite(full_ivw) else 0
    return max_delta, sign_flip


def steiger_proxy(beta_exp: np.ndarray, se_exp: np.ndarray, beta_out: np.ndarray, se_out: np.ndarray) -> tuple[bool, str]:
    r2_exp = np.square(beta_exp) / (np.square(beta_exp) + np.square(se_exp))
    r2_out = np.square(beta_out) / (np.square(beta_out) + np.square(se_out))
    support = bool(np.nanmean(r2_exp) > np.nanmean(r2_out))
    note = "proxy_without_eaf_or_sample_size"
    return support, note


def evaluate_direction(exposure_sig: pd.DataFrame, outcome: pd.DataFrame, direction: str, exposure_trait: str, outcome_trait: str, kb: int) -> dict:
    merged = prepare_merged(exposure_sig, outcome, kb)
    n_inst = int(len(merged))
    if n_inst < 3:
        return {
            "direction": direction,
            "exposure_trait": exposure_trait,
            "outcome_trait": outcome_trait,
            "n_instruments": n_inst,
            "cochrans_q": np.nan,
            "cochrans_q_df": max(n_inst - 1, 0),
            "cochrans_q_p": np.nan,
            "loo_max_abs_delta_ivw_beta": np.nan,
            "loo_sign_flip_count": 0,
            "steiger_proxy_supported": False,
            "steiger_proxy_note": "insufficient_instruments",
            "status": "insufficient_overlap",
        }

    beta_exp = merged["beta_exp"].to_numpy()
    se_exp = merged["se_exp"].to_numpy()
    beta_out = merged["beta_out"].to_numpy()
    se_out = merged["se_out"].to_numpy()

    ivw_beta, _, _ = ivw_mr(beta_exp, se_exp, beta_out, se_out)
    q, q_df, q_p = cochran_q(beta_exp, se_exp, beta_out, se_out)
    loo_max_delta, loo_flip = leave_one_out_ivw(beta_exp, se_exp, beta_out, se_out, ivw_beta)
    steiger_ok, steiger_note = steiger_proxy(beta_exp, se_exp, beta_out, se_out)
    return {
        "direction": direction,
        "exposure_trait": exposure_trait,
        "outcome_trait": outcome_trait,
        "n_instruments": n_inst,
        "cochrans_q": float(q),
        "cochrans_q_df": int(q_df),
        "cochrans_q_p": float(q_p) if np.isfinite(q_p) else np.nan,
        "loo_max_abs_delta_ivw_beta": float(loo_max_delta) if np.isfinite(loo_max_delta) else np.nan,
        "loo_sign_flip_count": int(loo_flip),
        "steiger_proxy_supported": bool(steiger_ok),
        "steiger_proxy_note": str(steiger_note),
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
    ibs_id = str(cfg["gwas"]["ibs_trait_id"])
    p_thr = float(cfg["analysis"]["mr_instrument_p"])
    kb = int(cfg["analysis"]["mr_clump_kb"])
    progress_every = int(cfg.get("runtime", {}).get("progress_every_row_groups", 20))

    rows: list[dict] = []
    psych_traits = [t for t in traits_cfg if t.get("group") == "psychiatric"]
    ibs_sig = load_trait_instruments(args.in_parquet, ibs_id, p_thr, progress_every)

    for trait in psych_traits:
        t_id = str(trait["trait_id"])
        print(f"[mr_sensitivity] trait={t_id} stage=load_instruments", flush=True)
        psych_sig = load_trait_instruments(args.in_parquet, t_id, p_thr, progress_every)
        ibs_for_psych = load_outcome_for_positions(args.in_parquet, ibs_id, positions_by_chr(psych_sig), progress_every)
        psych_for_ibs = load_outcome_for_positions(args.in_parquet, t_id, positions_by_chr(ibs_sig), progress_every)

        r1 = evaluate_direction(psych_sig, ibs_for_psych, f"{t_id}->IBS", t_id, ibs_id, kb)
        rows.append(r1)
        r2 = evaluate_direction(ibs_sig, psych_for_ibs, f"IBS->{t_id}", ibs_id, t_id, kb)
        rows.append(r2)

    out = pd.DataFrame(rows)
    out.to_csv(args.out_table, sep="\t", index=False)


if __name__ == "__main__":
    main()
