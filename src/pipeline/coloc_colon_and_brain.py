from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import yaml

from src.utils.stats import coloc_posteriors, wakefield_abf


OUTPUT_COLUMNS = [
    "locus_id",
    "lead_snp",
    "tissue",
    "gene_id",
    "n_variants",
    "top_eqtl_variant",
    "min_p_eqtl",
    "PP0",
    "PP1",
    "PP2",
    "PP3",
    "PP4",
    "mode",
    "exploratory",
]


def empty_output() -> pd.DataFrame:
    return pd.DataFrame(columns=OUTPUT_COLUMNS)


def normalize_leads(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty or not {"chr", "pos"}.issubset(frame.columns):
        return pd.DataFrame(columns=["lead_snp", "chr", "pos", "beta", "se", "p"])
    out = frame.copy()
    out["chr"] = out["chr"].astype(str)
    out["chr"] = out["chr"].where(out["chr"].str.startswith("chr"), "chr" + out["chr"])
    out["pos"] = pd.to_numeric(out["pos"], errors="coerce")
    out["beta"] = pd.to_numeric(out.get("beta", np.nan), errors="coerce")
    out["se"] = pd.to_numeric(out.get("se", np.nan), errors="coerce")
    out["p"] = pd.to_numeric(out.get("p", np.nan), errors="coerce")
    if "lead_snp" not in out.columns:
        rs = out.get("rsid", pd.Series([""] * len(out))).astype(str)
        out["lead_snp"] = rs.replace({"": np.nan}).fillna(out["chr"] + ":" + out["pos"].astype("Int64").astype(str))
    out = out.dropna(subset=["chr", "pos"]).copy()
    out["pos"] = out["pos"].astype(int)
    out["lead_snp"] = out["lead_snp"].astype(str)
    out = out.sort_values("p", na_position="last").drop_duplicates(subset=["chr", "pos"], keep="first")
    return out[["lead_snp", "chr", "pos", "beta", "se", "p"]]


def load_leads_from_gwas(gwas_path: Path, ibs_trait_id: str, progress_every: int) -> pd.DataFrame:
    pf = pq.ParquetFile(gwas_path)
    best_by_chr: dict[str, dict] = {}
    for rg in range(pf.metadata.num_row_groups):
        chunk = pf.read_row_group(rg, columns=["trait_id", "chr", "pos", "beta", "se", "p", "rsid"]).to_pandas()
        chunk = chunk.loc[chunk["trait_id"] == ibs_trait_id].copy()
        if chunk.empty:
            continue
        chunk["pos"] = pd.to_numeric(chunk["pos"], errors="coerce")
        chunk["p"] = pd.to_numeric(chunk["p"], errors="coerce")
        chunk = chunk.dropna(subset=["chr", "pos", "p"])
        if chunk.empty:
            continue
        chunk["chr"] = chunk["chr"].astype(str)
        chunk["chr"] = chunk["chr"].where(chunk["chr"].str.startswith("chr"), "chr" + chunk["chr"])
        chunk["beta"] = pd.to_numeric(chunk["beta"], errors="coerce")
        chunk["se"] = pd.to_numeric(chunk["se"], errors="coerce")
        for row in chunk.itertuples(index=False):
            chrom = str(row.chr)
            pvalue = float(row.p)
            prev = best_by_chr.get(chrom)
            if prev is None or pvalue < prev["p"]:
                rsid = str(getattr(row, "rsid", ""))
                if rsid and rsid != "nan":
                    lead_snp = rsid
                else:
                    lead_snp = f"{chrom}:{int(row.pos)}"
                best_by_chr[chrom] = {
                    "lead_snp": lead_snp,
                    "chr": chrom,
                    "pos": int(row.pos),
                    "beta": float(row.beta) if np.isfinite(row.beta) else np.nan,
                    "se": float(row.se) if np.isfinite(row.se) else np.nan,
                    "p": pvalue,
                }
        if (rg + 1) % progress_every == 0:
            print(
                f"[coloc] load_leads row_groups={rg + 1}/{pf.metadata.num_row_groups} best_chr={len(best_by_chr)}",
                flush=True,
            )
    if not best_by_chr:
        return pd.DataFrame(columns=["lead_snp", "chr", "pos", "beta", "se", "p"])
    out = pd.DataFrame(best_by_chr.values()).sort_values("p").head(50)
    return out[["lead_snp", "chr", "pos", "beta", "se", "p"]]


def load_leads(path: Path, gwas_path: Path, ibs_trait_id: str, progress_every: int) -> pd.DataFrame:
    if path.exists():
        try:
            leads = pd.read_csv(path, sep="\t")
            out = normalize_leads(leads)
            if not out.empty:
                return out
        except Exception:
            pass
    return load_leads_from_gwas(gwas_path, ibs_trait_id, progress_every)


def scan_brain_tissues(eqtl_path: Path) -> list[str]:
    pf = pq.ParquetFile(eqtl_path)
    tissues: set[str] = set()
    for rg in range(pf.metadata.num_row_groups):
        chunk = pf.read_row_group(rg, columns=["tissue"]).to_pandas()
        vals = chunk["tissue"].dropna().astype(str).unique().tolist()
        tissues.update(v for v in vals if "brain" in v.lower())
    return sorted(tissues)


def build_windows(leads: pd.DataFrame, window_bp: int) -> dict[str, list[dict]]:
    windows: dict[str, list[dict]] = {}
    for row in leads.itertuples(index=False):
        chrom = str(row.chr)
        pos = int(row.pos)
        windows.setdefault(chrom, []).append(
            {
                "locus_id": f"{chrom}:{pos}",
                "lead_snp": str(row.lead_snp),
                "start": pos - window_bp,
                "end": pos + window_bp,
                "lead_beta": float(row.beta) if np.isfinite(row.beta) else np.nan,
                "lead_se": float(row.se) if np.isfinite(row.se) else np.nan,
            }
        )
    return windows


def collect_eqtl_in_windows(eqtl_path: Path, tissues: set[str], windows: dict[str, list[dict]], progress_every: int) -> pd.DataFrame:
    cols = ["chr", "pos", "variant", "gene_id", "tissue", "beta", "se", "pvalue"]
    pf = pq.ParquetFile(eqtl_path)
    out_chunks: list[pd.DataFrame] = []
    valid_chr = set(windows.keys())
    for rg in range(pf.metadata.num_row_groups):
        chunk = pf.read_row_group(rg, columns=cols).to_pandas()
        chunk = chunk.loc[chunk["tissue"].isin(tissues) & chunk["chr"].isin(valid_chr)].copy()
        if chunk.empty:
            continue
        chunk["pos"] = pd.to_numeric(chunk["pos"], errors="coerce")
        chunk["beta"] = pd.to_numeric(chunk["beta"], errors="coerce")
        chunk["se"] = pd.to_numeric(chunk["se"], errors="coerce")
        chunk["pvalue"] = pd.to_numeric(chunk["pvalue"], errors="coerce")
        chunk = chunk.dropna(subset=["chr", "pos", "variant", "gene_id", "tissue", "beta", "se", "pvalue"])
        if chunk.empty:
            continue
        chunk["pos"] = chunk["pos"].astype(int)
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
            print(
                f"[coloc] eqtl_scan row_groups={rg + 1}/{pf.metadata.num_row_groups} collected_chunks={len(out_chunks)}",
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


def build_keys_by_chr(eqtl_local: pd.DataFrame) -> dict[str, set[int]]:
    keys: dict[str, set[int]] = {}
    if eqtl_local.empty:
        return keys
    for chrom, sub in eqtl_local.groupby("chr"):
        keys[str(chrom)] = set(sub["pos"].astype(int).tolist())
    return keys


def collect_gwas_map(
    gwas_path: Path,
    ibs_trait_id: str,
    keys_by_chr: dict[str, set[int]],
    progress_every: int,
) -> dict[tuple[str, int], tuple[float, float, float]]:
    if not keys_by_chr:
        return {}
    pf = pq.ParquetFile(gwas_path)
    best: dict[tuple[str, int], tuple[float, float, float]] = {}
    cols = ["trait_id", "chr", "pos", "beta", "se", "p"]
    valid_chr = set(keys_by_chr.keys())
    for rg in range(pf.metadata.num_row_groups):
        chunk = pf.read_row_group(rg, columns=cols).to_pandas()
        chunk = chunk.loc[(chunk["trait_id"] == ibs_trait_id) & chunk["chr"].isin(valid_chr)].copy()
        if chunk.empty:
            continue
        chunk["pos"] = pd.to_numeric(chunk["pos"], errors="coerce")
        chunk["beta"] = pd.to_numeric(chunk["beta"], errors="coerce")
        chunk["se"] = pd.to_numeric(chunk["se"], errors="coerce")
        chunk["p"] = pd.to_numeric(chunk["p"], errors="coerce")
        chunk = chunk.dropna(subset=["chr", "pos", "beta", "se", "p"])
        if chunk.empty:
            continue
        chunk["pos"] = chunk["pos"].astype(int)
        for chrom, sub in chunk.groupby("chr"):
            pos_set = keys_by_chr.get(str(chrom), set())
            sub = sub.loc[sub["pos"].isin(pos_set)]
            if sub.empty:
                continue
            for row in sub.itertuples(index=False):
                key = (str(row.chr), int(row.pos))
                rec = best.get(key)
                pvalue = float(row.p)
                if rec is None or pvalue < rec[2]:
                    best[key] = (float(row.beta), float(row.se), pvalue)
        if (rg + 1) % progress_every == 0:
            print(
                f"[coloc] gwas_scan row_groups={rg + 1}/{pf.metadata.num_row_groups} matched_sites={len(best)}",
                flush=True,
            )
    return best


def compute_coloc_for_group(
    sub: pd.DataFrame,
    gwas_map: dict[tuple[str, int], tuple[float, float, float]],
    lead_beta: float,
    lead_se: float,
    w: float,
    p1: float,
    p2: float,
    p12: float,
) -> tuple[dict, str]:
    abf_eqtl = wakefield_abf(sub["beta"].to_numpy(), sub["se"].to_numpy(), w)
    beta_gwas = np.full(len(sub), np.nan, dtype=float)
    se_gwas = np.full(len(sub), np.nan, dtype=float)
    for i, (chrom, pos) in enumerate(zip(sub["chr"].astype(str), sub["pos"].astype(int))):
        rec = gwas_map.get((chrom, pos))
        if rec is None:
            continue
        beta_gwas[i] = rec[0]
        se_gwas[i] = rec[1]
    have_variant_gwas = int(np.isfinite(beta_gwas).sum())

    if have_variant_gwas >= 3:
        abf_gwas = np.ones(len(sub), dtype=float)
        mask = np.isfinite(beta_gwas) & np.isfinite(se_gwas) & (se_gwas > 0)
        if mask.any():
            abf_gwas[mask] = wakefield_abf(beta_gwas[mask], se_gwas[mask], w)
        mode = "full"
    else:
        if np.isfinite(lead_beta) and np.isfinite(lead_se) and lead_se > 0:
            const_abf = float(wakefield_abf(np.array([lead_beta]), np.array([lead_se]), w)[0])
        else:
            const_abf = 1.0
        abf_gwas = np.full(len(sub), const_abf, dtype=float)
        mode = "proxy_lead_only"

    pp = coloc_posteriors(abf_eqtl, abf_gwas, p1, p2, p12)
    top_idx = int(np.nanargmax(abf_eqtl)) if len(abf_eqtl) else 0
    top_variant = str(sub.iloc[top_idx]["variant"]) if len(sub) else ""

    out = {
        "n_variants": int(len(sub)),
        "top_eqtl_variant": top_variant,
        "min_p_eqtl": float(sub["pvalue"].min()) if len(sub) else np.nan,
        "PP0": pp["PP0"],
        "PP1": pp["PP1"],
        "PP2": pp["PP2"],
        "PP3": pp["PP3"],
        "PP4": pp["PP4"],
    }
    return out, mode


def build_table(eqtl_local: pd.DataFrame, gwas_map: dict[tuple[str, int], tuple[float, float, float]], tissues: list[str], cfg: dict) -> pd.DataFrame:
    if eqtl_local.empty:
        return empty_output()
    subset = eqtl_local.loc[eqtl_local["tissue"].isin(tissues)].copy()
    if subset.empty:
        return empty_output()

    w = float(cfg["analysis"]["wakefield_W"])
    priors = cfg["analysis"]["coloc_priors"]
    p1 = float(priors["p1"])
    p2 = float(priors["p2"])
    p12 = float(priors["p12"])

    rows: list[dict] = []
    grouped = subset.groupby(["locus_id", "lead_snp", "tissue", "gene_id"], sort=False)
    for (locus_id, lead_snp, tissue, gene_id), sub in grouped:
        if len(sub) < 2:
            continue
        lead_beta = float(sub["lead_beta"].iloc[0]) if "lead_beta" in sub.columns else np.nan
        lead_se = float(sub["lead_se"].iloc[0]) if "lead_se" in sub.columns else np.nan
        metrics, mode = compute_coloc_for_group(sub, gwas_map, lead_beta, lead_se, w, p1, p2, p12)
        rows.append(
            {
                "locus_id": str(locus_id),
                "lead_snp": str(lead_snp),
                "tissue": str(tissue),
                "gene_id": str(gene_id),
                **metrics,
                "mode": mode,
                "exploratory": mode != "full",
            }
        )
    if not rows:
        return empty_output()
    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--eqtl", required=True)
    parser.add_argument("--gwas", required=True)
    parser.add_argument("--out-colon", required=True)
    parser.add_argument("--out-brain", required=True)
    args = parser.parse_args()

    cfg = yaml.safe_load(open(args.config, "r", encoding="utf-8"))
    ibs_trait_id = cfg["gwas"]["ibs_trait_id"]
    progress_every = int(cfg.get("runtime", {}).get("progress_every_row_groups", 20))
    leads = load_leads(Path("data_raw/ibsgwas_leads.tsv"), Path(args.gwas), ibs_trait_id, progress_every)
    if leads.empty:
        empty_output().to_csv(args.out_colon, sep="\t", index=False)
        empty_output().to_csv(args.out_brain, sep="\t", index=False)
        return

    window_bp = int(cfg["analysis"]["window_bp"])
    windows = build_windows(leads, window_bp)
    colon_tissues = ["colon_sigmoid", "colon_transverse"]
    brain_tissues = scan_brain_tissues(Path(args.eqtl))
    tissues_all = set(colon_tissues + brain_tissues)

    eqtl_local = collect_eqtl_in_windows(Path(args.eqtl), tissues_all, windows, progress_every)
    keys_by_chr = build_keys_by_chr(eqtl_local)
    gwas_map = collect_gwas_map(Path(args.gwas), ibs_trait_id, keys_by_chr, progress_every)

    coloc_colon = build_table(eqtl_local, gwas_map, colon_tissues, cfg)
    coloc_brain = build_table(eqtl_local, gwas_map, brain_tissues, cfg)
    coloc_colon.to_csv(args.out_colon, sep="\t", index=False)
    coloc_brain.to_csv(args.out_brain, sep="\t", index=False)


if __name__ == "__main__":
    main()
