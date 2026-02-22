from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def read_h5ad_with_anndata(path: Path):
    import anndata as ad

    return ad.read_h5ad(path, backed="r")


def choose_celltype_column(columns: list[str]) -> str | None:
    preferences = ["cell_type", "Cell_Type", "type", "lineage", "cluster", "leiden"]
    for col in preferences:
        if col in columns:
            return col
    return None


def get_vector(adata, gene: str, max_cells: int) -> tuple[np.ndarray, np.ndarray]:
    n = adata.n_obs
    idx = np.arange(min(n, max_cells))
    x = adata[idx, gene].X
    try:
        vec = x.toarray().ravel()
    except Exception:
        vec = np.asarray(x).ravel()
    return idx, vec


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--coloc-table", required=True)
    parser.add_argument("--out-table", required=True)
    args = parser.parse_args()

    cfg = yaml.safe_load(open(args.config, "r", encoding="utf-8"))

    coloc = pd.read_csv(args.coloc_table, sep="\t") if Path(args.coloc_table).exists() else pd.DataFrame()
    if not coloc.empty and "PP4" in coloc.columns:
        top = coloc.sort_values("PP4", ascending=False)
        genes = top.loc[top["PP4"] >= 0.5, "gene_id"].astype(str).unique().tolist()
        if not genes:
            genes = top["gene_id"].astype(str).head(20).tolist()
    else:
        genes = []

    rows: list[dict] = []
    max_cells = 10000

    for raw_path in cfg["sc"]["inputs"]:
        path = Path(raw_path)
        dataset_id = path.stem
        if not path.exists():
            rows.append({"dataset_id": dataset_id, "gene_id": "NA", "cell_type": "NA", "pct_expressing": np.nan, "mean_expr": np.nan, "status": "missing"})
            continue

        try:
            adata = read_h5ad_with_anndata(path)
        except Exception as exc:
            rows.append({"dataset_id": dataset_id, "gene_id": "NA", "cell_type": "NA", "pct_expressing": np.nan, "mean_expr": np.nan, "status": f"unreadable: {exc}"})
            continue

        var_names = set(map(str, adata.var_names))
        obs_cols = list(map(str, adata.obs.columns))
        cell_col = choose_celltype_column(obs_cols)
        if cell_col is None:
            rows.append({"dataset_id": dataset_id, "gene_id": "NA", "cell_type": "NA", "pct_expressing": np.nan, "mean_expr": np.nan, "status": "no_celltype_column"})
            adata.file.close()
            continue

        obs = adata.obs[[cell_col]].copy().iloc[:max_cells]
        obs[cell_col] = obs[cell_col].astype(str)
        groups = obs[cell_col].value_counts().head(20).index.tolist()

        selected_genes = [g for g in genes if g in var_names]
        if not selected_genes:
            rows.append({"dataset_id": dataset_id, "gene_id": "NA", "cell_type": "NA", "pct_expressing": np.nan, "mean_expr": np.nan, "status": "no_overlap_genes"})
            adata.file.close()
            continue

        for gene in selected_genes:
            idx, vec = get_vector(adata, gene, max_cells)
            local_obs = obs.iloc[: len(vec)].copy()
            local_obs["expr"] = vec
            for grp in groups:
                grp_vec = local_obs.loc[local_obs[cell_col] == grp, "expr"]
                if grp_vec.empty:
                    continue
                rows.append(
                    {
                        "dataset_id": dataset_id,
                        "gene_id": gene,
                        "cell_type": grp,
                        "pct_expressing": float((grp_vec > 0).mean()),
                        "mean_expr": float(grp_vec.mean()),
                        "status": "ok",
                    }
                )

        adata.file.close()

    out = pd.DataFrame(rows)
    if out.empty:
        out = pd.DataFrame(columns=["dataset_id", "gene_id", "cell_type", "pct_expressing", "mean_expr", "status"])
    out.to_csv(args.out_table, sep="\t", index=False)


if __name__ == "__main__":
    main()
