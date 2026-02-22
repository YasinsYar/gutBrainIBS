from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import yaml
from scipy.stats import pearsonr

from src.utils.stats import fdr_bh


def prune_best(best: dict[tuple[str, int], tuple[float, float, float]], max_rows: int) -> dict[tuple[str, int], tuple[float, float, float]]:
    if len(best) <= max_rows:
        return best
    keys = sorted(best, key=lambda k: best[k][2])[:max_rows]
    return {k: best[k] for k in keys}


def load_trait_pool(parquet_path: str, trait_id: str, p_ceiling: float, max_rows: int, progress_every: int) -> pd.DataFrame:
    pf = pq.ParquetFile(parquet_path)
    best: dict[tuple[str, int], tuple[float, float, float]] = {}
    for rg in range(pf.metadata.num_row_groups):
        chunk = pf.read_row_group(rg, columns=["trait_id", "chr", "pos", "beta", "se", "p"]).to_pandas()
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
        chunk = chunk.loc[chunk["p"] <= p_ceiling]
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
        if len(best) > max_rows * 2:
            best = prune_best(best, max_rows)
        if (rg + 1) % progress_every == 0:
            print(
                f"[genetic_correlation] trait={trait_id} row_groups={rg + 1}/{pf.metadata.num_row_groups} unique_sites={len(best)}",
                flush=True,
            )

    best = prune_best(best, max_rows)
    if not best:
        return pd.DataFrame(columns=["chr", "pos", "beta", "se", "p"])
    rows = [{"chr": k[0], "pos": k[1], "beta": v[0], "se": v[1], "p": v[2]} for k, v in best.items()]
    return pd.DataFrame(rows)


def merge_unique(ibs: pd.DataFrame, comp: pd.DataFrame) -> pd.DataFrame:
    if ibs.empty or comp.empty:
        return pd.DataFrame()
    ibs_idx = ibs.set_index(["chr", "pos"])
    comp_idx = comp.set_index(["chr", "pos"])
    common = ibs_idx.index.intersection(comp_idx.index)
    if len(common) == 0:
        return pd.DataFrame()
    ibs_take = ibs_idx.loc[common][["beta", "se"]].rename(columns={"beta": "beta_ibs", "se": "se_ibs"})
    comp_take = comp_idx.loc[common][["beta", "se"]].rename(columns={"beta": "beta_trait", "se": "se_trait"})
    out = ibs_take.join(comp_take, how="inner").reset_index()
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--in-parquet", required=True)
    parser.add_argument("--out-table", required=True)
    args = parser.parse_args()

    cfg = yaml.safe_load(open(args.config, "r", encoding="utf-8"))
    traits_cfg = yaml.safe_load(open(cfg["gwas"]["traits"], "r", encoding="utf-8"))["traits"]

    ibs_id = cfg["gwas"]["ibs_trait_id"]
    p_grid = sorted(float(v) for v in cfg["analysis"].get("rg_proxy_p_grid", [1.0e-4, 1.0e-3, 1.0e-2, 5.0e-2, 1.0]))
    max_rows = int(cfg["analysis"].get("rg_proxy_max_rows", 2_000_000))
    min_overlap = int(cfg["analysis"].get("rg_proxy_min_overlap", 30))
    p_ceiling = float(max(p_grid))
    progress_every = int(cfg.get("runtime", {}).get("progress_every_row_groups", 20))

    psych_traits = [t for t in traits_cfg if t.get("group") == "psychiatric"]
    rows = []

    ibs_pool = load_trait_pool(args.in_parquet, ibs_id, p_ceiling, max_rows, progress_every)

    for trait in psych_traits:
        t_id = trait["trait_id"]
        comp_pool = load_trait_pool(args.in_parquet, t_id, p_ceiling, max_rows, progress_every)
        chosen = None
        chosen_thr = None
        for thr in p_grid:
            ibs = ibs_pool.loc[ibs_pool["p"] <= thr, ["chr", "pos", "beta", "se", "p"]]
            comp = comp_pool.loc[comp_pool["p"] <= thr, ["chr", "pos", "beta", "se", "p"]]
            merged = merge_unique(ibs, comp)
            if len(merged) >= min_overlap:
                chosen = merged
                chosen_thr = thr
                break

        if chosen is None:
            rows.append(
                {
                    "trait_id": t_id,
                    "trait_label": trait["trait_label"],
                    "n_overlap": 0,
                    "rg_proxy": np.nan,
                    "p_proxy": np.nan,
                    "method": "proxy_z_correlation",
                    "p_threshold_used": np.nan,
                }
            )
            continue

        z1 = chosen["beta_ibs"] / chosen["se_ibs"]
        z2 = chosen["beta_trait"] / chosen["se_trait"]
        r, p = pearsonr(z1, z2)
        rows.append(
            {
                "trait_id": t_id,
                "trait_label": trait["trait_label"],
                "n_overlap": int(len(chosen)),
                "rg_proxy": float(r),
                "p_proxy": float(p),
                "method": "proxy_z_correlation",
                "p_threshold_used": float(chosen_thr),
            }
        )

    out = pd.DataFrame(rows)
    mask = out["p_proxy"].notna()
    out["q_fdr"] = np.nan
    if mask.any():
        out.loc[mask, "q_fdr"] = fdr_bh(out.loc[mask, "p_proxy"].to_numpy())
    out.to_csv(args.out_table, sep="\t", index=False)


if __name__ == "__main__":
    main()
