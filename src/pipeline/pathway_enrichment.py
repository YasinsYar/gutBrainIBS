from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import hypergeom

from src.utils.stats import fdr_bh


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--coloc-colon", required=True)
    parser.add_argument("--coloc-brain", required=True)
    parser.add_argument("--mr", required=True)
    parser.add_argument("--out-table", required=True)
    args = parser.parse_args()

    cfg = yaml.safe_load(open(args.config, "r", encoding="utf-8"))

    coloc_brain = pd.read_csv(args.coloc_brain, sep="\t") if Path(args.coloc_brain).exists() else pd.DataFrame()
    selected = set(coloc_brain.loc[coloc_brain.get("PP4", pd.Series(dtype=float)) >= 0.5, "gene_id"].astype(str))

    pathway_file = cfg.get("pathways", {}).get("gene_sets_tsv", "")
    rows: list[dict] = []

    if not pathway_file or not Path(pathway_file).exists():
        rows.append(
            {
                "pathway_id": "NA",
                "description": "No pathway database provided",
                "n_input": len(selected),
                "overlap": 0,
                "p_value": np.nan,
                "q_value": np.nan,
                "method": "hypergeometric",
                "status": "not_run_missing_gene_sets",
            }
        )
        pd.DataFrame(rows).to_csv(args.out_table, sep="\t", index=False)
        return

    genesets = pd.read_csv(pathway_file, sep="\t")
    if not {"pathway_id", "gene_id"}.issubset(genesets.columns):
        rows.append(
            {
                "pathway_id": "NA",
                "description": "gene_sets_tsv schema mismatch",
                "n_input": len(selected),
                "overlap": 0,
                "p_value": np.nan,
                "q_value": np.nan,
                "method": "hypergeometric",
                "status": "schema_mismatch",
            }
        )
        pd.DataFrame(rows).to_csv(args.out_table, sep="\t", index=False)
        return

    universe = set(genesets["gene_id"].astype(str).unique())
    hits = set(g for g in selected if g in universe)
    M = len(universe)
    n = len(hits)

    pvals = []
    temp = []
    for pid, sub in genesets.groupby("pathway_id"):
        genes = set(sub["gene_id"].astype(str))
        k = len(hits.intersection(genes))
        N = len(genes)
        if n == 0 or M == 0 or N == 0:
            p = np.nan
        else:
            p = hypergeom.sf(k - 1, M, N, n)
        pvals.append(p)
        temp.append(
            {
                "pathway_id": pid,
                "description": str(sub.get("description", pd.Series([pid])).iloc[0]),
                "n_input": n,
                "overlap": k,
                "p_value": p,
                "method": "hypergeometric",
                "status": "ok",
            }
        )

    out = pd.DataFrame(temp)
    mask = out["p_value"].notna()
    out["q_value"] = np.nan
    if mask.any():
        out.loc[mask, "q_value"] = fdr_bh(out.loc[mask, "p_value"].to_numpy())

    out = out.sort_values(["q_value", "p_value"], na_position="last")
    out.to_csv(args.out_table, sep="\t", index=False)


if __name__ == "__main__":
    main()
