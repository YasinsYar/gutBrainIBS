from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def read_tsv(path: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(p, sep="\t")
    except Exception:
        return pd.DataFrame()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coloc-summary", required=True)
    parser.add_argument("--mr-robustness", required=True)
    parser.add_argument("--harmonization-compare", required=True)
    parser.add_argument("--out-table", required=True)
    args = parser.parse_args()

    coloc = read_tsv(args.coloc_summary)
    mr = read_tsv(args.mr_robustness)
    harm = read_tsv(args.harmonization_compare)

    rows: list[dict] = []

    if coloc.empty:
        rows.append({"module": "coloc_sensitivity", "status": "missing_or_empty"})
    else:
        base = coloc.loc[coloc["scenario_id"].astype(str).str.endswith("_all")].copy()
        if base.empty:
            base = coloc.copy()
        rows.append(
            {
                "module": "coloc_sensitivity",
                "status": "ok",
                "n_scenarios": int(len(coloc)),
                "max_pp4": float(coloc["pp4_max"].max()) if "pp4_max" in coloc.columns else None,
                "max_n_pp4_ge_0_7": int(coloc.get("n_pp4_ge_0_7", pd.Series(dtype=float)).max()) if "n_pp4_ge_0_7" in coloc.columns else 0,
                "top_changed_in_any_scenario": bool(coloc.get("top_hit_changed_vs_baseline", pd.Series(dtype=bool)).astype(bool).any())
                if "top_hit_changed_vs_baseline" in coloc.columns
                else False,
            }
        )

    if mr.empty:
        rows.append({"module": "mr_robustness", "status": "missing_or_empty"})
    else:
        ok = mr.loc[mr.get("status", pd.Series(dtype=str)) == "ok"].copy()
        rows.append(
            {
                "module": "mr_robustness",
                "status": "ok",
                "n_rows": int(len(mr)),
                "n_ok_rows": int(len(ok)),
                "n_q_p_lt_0_05": int((ok.get("cochrans_q_p", pd.Series(dtype=float)) < 0.05).sum()),
                "n_loo_sign_flip": int((ok.get("loo_sign_flip_count", pd.Series(dtype=float)) > 0).sum()),
                "n_steiger_proxy_supported": int(ok.get("steiger_proxy_supported", pd.Series(dtype=bool)).astype(bool).sum()),
            }
        )

    if harm.empty:
        rows.append({"module": "harmonization_strict_vs_nonstrict", "status": "missing_or_empty"})
    else:
        r = harm.iloc[0]
        rows.append(
            {
                "module": "harmonization_strict_vs_nonstrict",
                "status": str(r.get("status", "ok")),
                "n_added_pal_rows": int(r.get("n_added_palindromic_window_rows", 0)),
                "top_hit_changed": bool(r.get("top_hit_changed", False)),
                "strict_n_pp4_ge_0_7": int(r.get("strict_n_pp4_ge_0_7", 0)),
                "nonstrict_n_pp4_ge_0_7": int(r.get("nonstrict_n_pp4_ge_0_7", 0)),
                "strict_label": str(r.get("strict_label", "NA")),
                "nonstrict_label": str(r.get("nonstrict_label", "NA")),
            }
        )

    out = pd.DataFrame(rows)
    Path(args.out_table).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_table, sep="\t", index=False)


if __name__ == "__main__":
    main()
