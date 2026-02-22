from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import yaml


def to_bool(x) -> bool:
    if pd.isna(x):
        return False
    return bool(x)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--validation-json", required=True)
    parser.add_argument("--shared-metrics", required=True)
    parser.add_argument("--rg", required=True)
    parser.add_argument("--mr", required=True)
    parser.add_argument("--coloc-brain", required=True)
    parser.add_argument("--out-scorecard", required=True)
    parser.add_argument("--out-evidence", required=True)
    parser.add_argument("--out-qc-html", required=True)
    args = parser.parse_args()

    cfg = yaml.safe_load(open(args.config, "r", encoding="utf-8"))
    alpha = float(cfg["analysis"]["fdr_alpha"])

    validation = json.load(open(args.validation_json, "r", encoding="utf-8"))
    shared_metrics = json.load(open(args.shared_metrics, "r", encoding="utf-8")) if Path(args.shared_metrics).exists() else {}

    rg = pd.read_csv(args.rg, sep="\t") if Path(args.rg).exists() else pd.DataFrame()
    mr = pd.read_csv(args.mr, sep="\t") if Path(args.mr).exists() else pd.DataFrame()
    coloc = pd.read_csv(args.coloc_brain, sep="\t") if Path(args.coloc_brain).exists() else pd.DataFrame()

    crit_rg = bool((rg.get("q_fdr", pd.Series(dtype=float)) < alpha).any()) if not rg.empty else False
    crit_coloc = bool((coloc.get("PP4", pd.Series(dtype=float)) >= 0.7).any()) if not coloc.empty else False

    mr_ok = False
    if not mr.empty:
        rows = mr.dropna(subset=["ivw_p", "consistent_direction"]).copy()
        if not rows.empty:
            mr_ok = bool(((rows["ivw_p"] < alpha) & (rows["consistent_direction"].astype(bool)) & (~rows["pleiotropy_flag"].astype(bool))).any())

    if crit_rg and crit_coloc and mr_ok:
        label = "Supported gut-brain genetic contribution"
    elif crit_rg or crit_coloc or mr_ok:
        label = "Partial support"
    else:
        label = "Insufficient support"

    exploratory_coloc = False
    if not coloc.empty and "exploratory" in coloc.columns:
        exploratory_coloc = bool(coloc["exploratory"].astype(bool).all())

    scorecard = {
        "research_frame": "Biopsychosocial DGBI",
        "interpretation_label": label,
        "guardrail": "No claim of psychosomatic origin proven.",
        "criteria": {
            "genetic_correlation_supported": crit_rg,
            "strong_brain_coloc_supported": crit_coloc,
            "mr_directional_support": mr_ok,
        },
        "exploratory_coloc_only": exploratory_coloc,
        "validation": validation,
        "shared_specific_metrics": shared_metrics,
    }

    Path(args.out_scorecard).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_scorecard, "w", encoding="utf-8") as handle:
        json.dump(scorecard, handle, indent=2)

    evidence = pd.DataFrame(
        [
            {
                "criterion": "IBS-psychiatric genetic correlation",
                "met": crit_rg,
                "details": "Any psychiatric trait q_fdr < alpha",
            },
            {
                "criterion": "Strong brain colocalization",
                "met": crit_coloc,
                "details": "Any brain PP4 >= 0.7",
            },
            {
                "criterion": "Bidirectional MR support",
                "met": mr_ok,
                "details": "Any direction with IVW p < alpha, consistent direction, no pleiotropy flag",
            },
            {
                "criterion": "Final interpretation",
                "met": True,
                "details": label,
            },
        ]
    )
    evidence.to_csv(args.out_evidence, sep="\t", index=False)

    html = f"""
    <html>
      <head><title>Run QC</title></head>
      <body>
        <h1>Gut-Brain IBS Pipeline QC</h1>
        <p><b>Interpretation:</b> {label}</p>
        <p><b>Guardrail:</b> No claim of psychosomatic origin proven.</p>
        <h2>Input validation</h2>
        <pre>{json.dumps(validation, indent=2)}</pre>
        <h2>Shared-specific metrics</h2>
        <pre>{json.dumps(shared_metrics, indent=2)}</pre>
        <h2>Criteria</h2>
        <ul>
          <li>Genetic correlation support: {crit_rg}</li>
          <li>Strong brain coloc support: {crit_coloc}</li>
          <li>MR support: {mr_ok}</li>
          <li>Exploratory coloc only: {exploratory_coloc}</li>
        </ul>
      </body>
    </html>
    """
    Path(args.out_qc_html).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_qc_html).write_text(html, encoding="utf-8")


if __name__ == "__main__":
    main()
