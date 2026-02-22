from __future__ import annotations

import argparse

import yaml

from src.utils.duckdb_utils import connect_duckdb, sql_path


def heterogeneity_query(in_path: str) -> str:
    return f"""
    WITH dedup AS (
      SELECT
        CAST(variant AS VARCHAR) AS variant,
        CAST(gene_id AS VARCHAR) AS gene_id,
        CAST(tissue AS VARCHAR) AS tissue,
        TRY_CAST(beta AS DOUBLE) AS beta,
        TRY_CAST(se AS DOUBLE) AS se,
        TRY_CAST(pvalue AS DOUBLE) AS pvalue,
        row_number() OVER (
          PARTITION BY CAST(variant AS VARCHAR), CAST(gene_id AS VARCHAR), CAST(tissue AS VARCHAR)
          ORDER BY TRY_CAST(pvalue AS DOUBLE) ASC NULLS LAST
        ) AS rn
      FROM read_parquet({sql_path(in_path)})
      WHERE tissue IN ('colon_sigmoid', 'colon_transverse')
    ),
    best AS (
      SELECT variant, gene_id, tissue, beta, se
      FROM dedup
      WHERE rn = 1
    ),
    agg AS (
      SELECT
        variant,
        gene_id,
        max(CASE WHEN tissue = 'colon_sigmoid' THEN beta END) AS beta_colon_sigmoid,
        max(CASE WHEN tissue = 'colon_transverse' THEN beta END) AS beta_colon_transverse,
        max(CASE WHEN tissue = 'colon_sigmoid' THEN se END) AS se_colon_sigmoid,
        max(CASE WHEN tissue = 'colon_transverse' THEN se END) AS se_colon_transverse
      FROM best
      GROUP BY variant, gene_id
    ),
    calc AS (
      SELECT
        variant,
        gene_id,
        (beta_colon_sigmoid - beta_colon_transverse)
          / sqrt((se_colon_sigmoid * se_colon_sigmoid) + (se_colon_transverse * se_colon_transverse)) AS z_het
      FROM agg
      WHERE beta_colon_sigmoid IS NOT NULL
        AND beta_colon_transverse IS NOT NULL
        AND se_colon_sigmoid IS NOT NULL
        AND se_colon_transverse IS NOT NULL
        AND ((se_colon_sigmoid * se_colon_sigmoid) + (se_colon_transverse * se_colon_transverse)) > 0
    ),
    pvals AS (
      SELECT
        variant,
        gene_id,
        z_het,
        least(
          1.0,
          greatest(
            0.0,
            2.0 * (
              exp(-(abs(z_het) * abs(z_het)) / 2.0) / sqrt(2.0 * pi())
            ) * (
              (0.319381530 * t) +
              (-0.356563782 * t * t) +
              (1.781477937 * t * t * t) +
              (-1.821255978 * t * t * t * t) +
              (1.330274429 * t * t * t * t * t)
            )
          )
        ) AS p_het
      FROM calc
      CROSS JOIN LATERAL (
        SELECT 1.0 / (1.0 + 0.2316419 * abs(z_het)) AS t
      ) tt
    ),
    ranked AS (
      SELECT
        *,
        row_number() OVER (ORDER BY p_het ASC NULLS LAST) AS r,
        count(*) OVER () AS n
      FROM pvals
    ),
    raw_q AS (
      SELECT
        *,
        least(1.0, (p_het * n) / r) AS q_raw
      FROM ranked
    ),
    adjusted AS (
      SELECT
        *,
        min(q_raw) OVER (
          ORDER BY p_het DESC NULLS LAST
          ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW
        ) AS q_het
      FROM raw_q
    )
    SELECT
      variant,
      gene_id,
      z_het,
      p_het,
      least(1.0, q_het) AS q_het
    FROM adjusted
    ORDER BY p_het ASC NULLS LAST
    """


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--in-parquet", required=True)
    parser.add_argument("--out-table", required=True)
    args = parser.parse_args()

    cfg = yaml.safe_load(open(args.config, "r", encoding="utf-8"))
    con = connect_duckdb(cfg)
    try:
        con.execute(
            f"""
            COPY (
              {heterogeneity_query(args.in_parquet)}
            ) TO {sql_path(args.out_table)} (HEADER TRUE, DELIMITER '\\t');
            """
        )
    finally:
        con.close()


if __name__ == "__main__":
    main()
