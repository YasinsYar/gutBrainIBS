from __future__ import annotations

import argparse

import yaml

from src.utils.duckdb_utils import connect_duckdb, sql_path
from src.utils.io import write_json


def final_cte(in_path: str) -> str:
    return f"""
    WITH dedup AS (
      SELECT
        CAST(variant AS VARCHAR) AS variant,
        CAST(gene_id AS VARCHAR) AS gene_id,
        CAST(tissue AS VARCHAR) AS tissue,
        TRY_CAST(beta AS DOUBLE) AS beta,
        TRY_CAST(pvalue AS DOUBLE) AS pvalue,
        row_number() OVER (
          PARTITION BY CAST(variant AS VARCHAR), CAST(gene_id AS VARCHAR), CAST(tissue AS VARCHAR)
          ORDER BY TRY_CAST(pvalue AS DOUBLE) ASC NULLS LAST
        ) AS rn
      FROM read_parquet({sql_path(in_path)})
      WHERE tissue IN ('colon_sigmoid', 'colon_transverse')
    ),
    best AS (
      SELECT variant, gene_id, tissue, beta
      FROM dedup
      WHERE rn = 1
    ),
    agg AS (
      SELECT
        variant,
        gene_id,
        max(CASE WHEN tissue = 'colon_sigmoid' THEN beta END) AS beta_sigmoid,
        max(CASE WHEN tissue = 'colon_transverse' THEN beta END) AS beta_transverse
      FROM best
      GROUP BY variant, gene_id
    )
    SELECT
      variant,
      gene_id,
      (beta_sigmoid IS NOT NULL) AS in_sigmoid,
      (beta_transverse IS NOT NULL) AS in_transverse,
      beta_sigmoid,
      beta_transverse,
      abs(beta_sigmoid - beta_transverse) AS delta_beta_abs,
      CASE
        WHEN beta_sigmoid IS NOT NULL AND beta_transverse IS NOT NULL THEN sign(beta_sigmoid) = sign(beta_transverse)
        ELSE NULL
      END AS sign_concordant,
      CASE
        WHEN beta_sigmoid IS NOT NULL AND beta_transverse IS NOT NULL THEN 'shared'
        WHEN beta_sigmoid IS NOT NULL AND beta_transverse IS NULL THEN 'sigmoid_specific'
        WHEN beta_sigmoid IS NULL AND beta_transverse IS NOT NULL THEN 'transverse_specific'
        ELSE 'shared'
      END AS category
    FROM agg
    """


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--in-parquet", required=True)
    parser.add_argument("--out-table", required=True)
    parser.add_argument("--out-metrics", required=True)
    args = parser.parse_args()

    cfg = yaml.safe_load(open(args.config, "r", encoding="utf-8"))
    con = connect_duckdb(cfg)
    cte = final_cte(args.in_parquet)
    try:
        con.execute(
            f"""
            COPY (
              {cte}
            ) TO {sql_path(args.out_table)} (HEADER TRUE, DELIMITER '\\t');
            """
        )

        row = con.execute(
            f"""
            SELECT
              count(*) AS n_pairs_total,
              sum(CASE WHEN category = 'shared' THEN 1 ELSE 0 END) AS n_shared,
              sum(CASE WHEN category <> 'shared' THEN 1 ELSE 0 END) AS n_tissue_specific,
              sum(CASE WHEN category = 'sigmoid_specific' THEN 1 ELSE 0 END) AS n_sigmoid_specific,
              sum(CASE WHEN category = 'transverse_specific' THEN 1 ELSE 0 END) AS n_transverse_specific,
              avg(CASE WHEN category = 'shared' THEN CAST(sign_concordant AS DOUBLE) END) AS sign_concordance,
              avg(CASE WHEN category = 'shared' THEN delta_beta_abs END) AS delta_beta_mean,
              median(CASE WHEN category = 'shared' THEN delta_beta_abs END) AS delta_beta_median
            FROM ({cte}) t;
            """
        ).fetchone()
    finally:
        con.close()

    metrics = {
        "n_pairs_total": int(row[0] or 0),
        "n_shared": int(row[1] or 0),
        "n_tissue_specific": int(row[2] or 0),
        "n_sigmoid_specific": int(row[3] or 0),
        "n_transverse_specific": int(row[4] or 0),
        "sign_concordance": float(row[5] or 0.0),
        "delta_beta_mean": float(row[6] or 0.0),
        "delta_beta_median": float(row[7] or 0.0),
    }
    write_json(metrics, args.out_metrics)


if __name__ == "__main__":
    main()
