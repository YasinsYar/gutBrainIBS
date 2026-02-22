from __future__ import annotations

import argparse
import contextlib
import io
import json
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text
from matplotlib_venn import venn2
from upsetplot import UpSet, from_indicators


def read_tsv(path: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(p, sep="\t")
    except Exception:
        return pd.DataFrame()


def save_placeholder(path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.text(0.5, 0.5, "No data available", ha="center", va="center", fontsize=13)
    ax.set_title(title)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def style() -> None:
    sns.set_theme(style="whitegrid", context="talk")
    plt.rcParams["figure.facecolor"] = "#f7f8fa"
    plt.rcParams["axes.facecolor"] = "#f7f8fa"
    plt.rcParams["axes.edgecolor"] = "#4b5563"
    plt.rcParams["axes.titleweight"] = "bold"
    plt.rcParams["axes.labelcolor"] = "#111827"
    plt.rcParams["xtick.color"] = "#111827"
    plt.rcParams["ytick.color"] = "#111827"
    plt.rcParams["font.size"] = 12


def fig_upset(shared: pd.DataFrame, out_path: Path) -> None:
    required = {"in_sigmoid", "in_transverse"}
    if shared.empty or not required.issubset(shared.columns):
        save_placeholder(out_path, "Colon eQTL overlap (UpSet)")
        return
    inds = shared.loc[:, ["in_sigmoid", "in_transverse"]].fillna(False).astype(bool)
    series = from_indicators(["in_sigmoid", "in_transverse"], inds)
    fig = plt.figure(figsize=(10, 6))
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
        warnings.filterwarnings("ignore", category=UserWarning)
        UpSet(series, subset_size="count", show_counts=True, sort_by="cardinality").plot(fig=fig)
    plt.suptitle("Colon eQTL overlap between sigmoid and transverse", y=1.02)
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def fig_venn(shared: pd.DataFrame, out_path: Path) -> None:
    required = {"in_sigmoid", "in_transverse"}
    if shared.empty or not required.issubset(shared.columns):
        save_placeholder(out_path, "Colon eQTL overlap (Venn)")
        return
    sigmoid = set(shared.loc[shared["in_sigmoid"].astype(bool), "variant"].astype(str))
    transverse = set(shared.loc[shared["in_transverse"].astype(bool), "variant"].astype(str))
    fig, ax = plt.subplots(figsize=(7, 6))
    venn2([sigmoid, transverse], set_labels=("Sigmoid", "Transverse"), ax=ax)
    ax.set_title("Colon eQTL variant overlap")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def fig_beta_hex(shared: pd.DataFrame, out_path: Path) -> None:
    required = {"beta_sigmoid", "beta_transverse"}
    if shared.empty or not required.issubset(shared.columns):
        save_placeholder(out_path, "Beta concordance (hexbin)")
        return
    df = shared.dropna(subset=["beta_sigmoid", "beta_transverse"]).copy()
    if df.empty:
        save_placeholder(out_path, "Beta concordance (hexbin)")
        return
    x = df["beta_sigmoid"].to_numpy()
    y = df["beta_transverse"].to_numpy()
    lim = float(np.nanmax(np.abs(np.concatenate([x, y]))))
    fig, ax = plt.subplots(figsize=(8, 7))
    hb = ax.hexbin(x, y, gridsize=45, cmap="mako", mincnt=1)
    ax.plot([-lim, lim], [-lim, lim], linestyle="--", color="#111827", linewidth=1.2)
    corr = float(np.corrcoef(x, y)[0, 1]) if len(df) > 1 else np.nan
    ax.set_title(f"Sigmoid vs transverse effect sizes (r={corr:.3f})")
    ax.set_xlabel("beta_sigmoid")
    ax.set_ylabel("beta_transverse")
    cb = fig.colorbar(hb, ax=ax)
    cb.set_label("Pair density")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def fig_heterogeneity(het: pd.DataFrame, out_path: Path) -> None:
    required = {"z_het", "p_het", "q_het", "gene_id"}
    if het.empty or not required.issubset(het.columns):
        save_placeholder(out_path, "Heterogeneity volcano")
        return
    df = het.copy()
    df["neglog10p"] = -np.log10(np.clip(df["p_het"].astype(float), 1e-300, 1.0))
    df["sig_fdr"] = df["q_het"].astype(float) < 0.05
    fig, ax = plt.subplots(figsize=(10, 7))
    sns.scatterplot(
        data=df,
        x="z_het",
        y="neglog10p",
        hue="sig_fdr",
        palette={True: "#c81e1e", False: "#6b7280"},
        s=30,
        alpha=0.7,
        linewidth=0,
        ax=ax,
        legend=False,
    )
    ax.axhline(-np.log10(0.05), linestyle="--", color="#111827", linewidth=1)
    ax.set_title("Colon tissue heterogeneity (volcano)")
    ax.set_xlabel("z_het")
    ax.set_ylabel("-log10(p_het)")

    lab = df.sort_values("q_het", ascending=True).head(12).copy()
    texts = []
    for row in lab.itertuples(index=False):
        texts.append(ax.text(float(row.z_het), float(row.neglog10p), str(row.gene_id), fontsize=9, color="#111827"))
    if texts:
        with contextlib.redirect_stdout(io.StringIO()):
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#374151", lw=0.6))
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def fig_coloc_heatmap(coloc_brain: pd.DataFrame, out_path: Path) -> None:
    required = {"gene_id", "tissue", "PP4"}
    if coloc_brain.empty or not required.issubset(coloc_brain.columns):
        save_placeholder(out_path, "Brain colocalization heatmap")
        return
    df = coloc_brain.copy()
    top_genes = (
        df.groupby("gene_id", as_index=False)["PP4"]
        .max()
        .sort_values("PP4", ascending=False)
        .head(16)["gene_id"]
        .astype(str)
        .tolist()
    )
    sub = df.loc[df["gene_id"].astype(str).isin(top_genes)].copy()
    if sub.empty:
        save_placeholder(out_path, "Brain colocalization heatmap")
        return
    pivot = sub.pivot_table(index="tissue", columns="gene_id", values="PP4", aggfunc="max", fill_value=0.0)
    fig, ax = plt.subplots(figsize=(14, 8))
    sns.heatmap(pivot, cmap="rocket_r", vmin=0, vmax=1, linewidths=0.4, linecolor="#e5e7eb", ax=ax)
    ax.set_title("Brain colocalization PP4 heatmap (top genes)")
    ax.set_xlabel("gene_id")
    ax.set_ylabel("brain tissue")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def fig_coloc_locus_bubble(coloc_brain: pd.DataFrame, out_path: Path) -> None:
    required = {"locus_id", "tissue", "PP4", "gene_id"}
    if coloc_brain.empty or not required.issubset(coloc_brain.columns):
        save_placeholder(out_path, "Brain colocalization bubble map")
        return
    df = coloc_brain.copy()
    top_loci = df.groupby("locus_id", as_index=False)["PP4"].max().sort_values("PP4", ascending=False).head(12)["locus_id"].tolist()
    sub = df.loc[df["locus_id"].isin(top_loci)].copy()
    if sub.empty:
        save_placeholder(out_path, "Brain colocalization bubble map")
        return
    sub["locus_rank"] = pd.Categorical(sub["locus_id"], categories=top_loci, ordered=True).codes
    fig, ax = plt.subplots(figsize=(12, 7))
    sc = ax.scatter(
        sub["locus_rank"],
        sub["tissue"].astype(str),
        s=60 + 600 * sub["PP4"].clip(0, 1),
        c=sub["PP4"],
        cmap="viridis",
        alpha=0.75,
        edgecolors="#111827",
        linewidths=0.3,
    )
    ax.set_xticks(range(len(top_loci)))
    ax.set_xticklabels(top_loci, rotation=45, ha="right")
    ax.set_title("Brain colocalization across top IBS loci")
    ax.set_xlabel("locus_id")
    ax.set_ylabel("tissue")
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("PP4")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def fig_coloc_sensitivity(coloc_sens: pd.DataFrame, out_path: Path) -> None:
    required = {"window_bp", "p12", "mhc_mode", "n_pp4_ge_0_7"}
    if coloc_sens.empty or not required.issubset(coloc_sens.columns):
        save_placeholder(out_path, "Colocalization sensitivity heatmap")
        return
    base = coloc_sens.loc[coloc_sens["mhc_mode"].astype(str) == "all"].copy()
    if base.empty:
        save_placeholder(out_path, "Colocalization sensitivity heatmap")
        return
    pivot = base.pivot_table(index="window_bp", columns="p12", values="n_pp4_ge_0_7", aggfunc="max", fill_value=0)
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(pivot, annot=True, fmt=".0f", cmap="crest", cbar_kws={"label": "N(PP4 >= 0.7)"}, ax=ax)
    ax.set_title("Brain colocalization sensitivity: window x p12")
    ax.set_xlabel("p12 prior")
    ax.set_ylabel("window_bp")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def fig_mhc_sensitivity(coloc_sens: pd.DataFrame, out_path: Path) -> None:
    required = {"scenario_id", "mhc_mode", "n_pp4_ge_0_7", "pp4_max"}
    if coloc_sens.empty or not required.issubset(coloc_sens.columns):
        save_placeholder(out_path, "MHC sensitivity")
        return
    baseline = coloc_sens.loc[coloc_sens["scenario_id"].astype(str).str.contains("w500000_p12_1e-05_")].copy()
    if baseline.empty:
        baseline = coloc_sens.loc[coloc_sens["mhc_mode"].astype(str).isin(["all", "exclude_mhc", "chr6_only"])].copy()
    if baseline.empty:
        save_placeholder(out_path, "MHC sensitivity")
        return
    baseline = baseline.sort_values("mhc_mode")
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.barplot(data=baseline, x="mhc_mode", y="n_pp4_ge_0_7", hue="mhc_mode", palette="Set2", dodge=False, legend=False, ax=ax)
    for i, row in enumerate(baseline.itertuples(index=False)):
        ax.text(i, float(row.n_pp4_ge_0_7) + 0.5, f"maxPP4={float(row.pp4_max):.3f}", ha="center", fontsize=10)
    ax.set_title("MHC sensitivity for strong brain colocalization")
    ax.set_xlabel("scenario")
    ax.set_ylabel("N(PP4 >= 0.7)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def fig_rg_panel(rg: pd.DataFrame, out_path: Path) -> None:
    required = {"trait_id", "rg_proxy", "q_fdr"}
    if rg.empty or not required.issubset(rg.columns):
        save_placeholder(out_path, "Genetic correlation panel")
        return
    df = rg.copy().sort_values("rg_proxy", ascending=False)
    df["sig"] = df["q_fdr"].astype(float) < 0.05
    fig, ax = plt.subplots(figsize=(8, 6))
    for i, row in enumerate(df.itertuples(index=False)):
        color = "#0f766e" if bool(row.sig) else "#9ca3af"
        ax.plot([0, float(row.rg_proxy)], [i, i], color=color, linewidth=3, alpha=0.8)
        ax.scatter([float(row.rg_proxy)], [i], color=color, s=120)
        ax.text(float(row.rg_proxy), i + 0.1, f"q={float(row.q_fdr):.3g}", fontsize=10)
    ax.axvline(0, color="#111827", linestyle="--", linewidth=1)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["trait_id"].astype(str).tolist())
    ax.set_xlabel("rg_proxy")
    ax.set_title("IBS vs psychiatric traits (proxy genetic correlation)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def fig_mr_robustness(mr: pd.DataFrame, mr_rob: pd.DataFrame, out_path: Path) -> None:
    required_mr = {"direction", "ivw_beta", "ivw_se", "ivw_p", "status"}
    if mr.empty or not required_mr.issubset(mr.columns):
        save_placeholder(out_path, "MR forest with robustness")
        return
    plot_df = mr.loc[mr["status"].astype(str) == "ok"].copy()
    if plot_df.empty:
        save_placeholder(out_path, "MR forest with robustness")
        return
    if not mr_rob.empty and {"direction", "cochrans_q_p", "loo_sign_flip_count"}.issubset(mr_rob.columns):
        plot_df = plot_df.merge(mr_rob[["direction", "cochrans_q_p", "loo_sign_flip_count"]], on="direction", how="left")
    else:
        plot_df["cochrans_q_p"] = np.nan
        plot_df["loo_sign_flip_count"] = 0
    plot_df = plot_df.sort_values("ivw_beta").reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(11, 7))
    y = np.arange(len(plot_df))
    is_sig = plot_df["ivw_p"].astype(float) < 0.05
    colors = np.where(is_sig, "#1d4ed8", "#9ca3af")
    x = plot_df["ivw_beta"].astype(float).to_numpy()
    xerr = 1.96 * plot_df["ivw_se"].astype(float).to_numpy()
    for i in range(len(plot_df)):
        ax.errorbar(x[i], y[i], xerr=xerr[i], fmt="o", color=colors[i], ecolor=colors[i], capsize=3, markersize=8)
        q = plot_df["cochrans_q_p"].iloc[i]
        loo = int(plot_df["loo_sign_flip_count"].iloc[i]) if not pd.isna(plot_df["loo_sign_flip_count"].iloc[i]) else 0
        qtxt = f"Q={q:.2g}" if not pd.isna(q) else "Q=NA"
        ax.text(x[i], y[i] + 0.15, f"{qtxt}; LOOflip={loo}", fontsize=9)
    ax.axvline(0, color="#111827", linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["direction"].astype(str).tolist())
    ax.set_xlabel("IVW beta (95% CI)")
    ax.set_title("Bidirectional MR with robustness annotations")
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def build_tables(
    coloc_brain: pd.DataFrame,
    mr: pd.DataFrame,
    mr_rob: pd.DataFrame,
    coloc_sens: pd.DataFrame,
    harm_sens: pd.DataFrame,
    table_dir: Path,
    html_path: Path,
    md_path: Path,
) -> None:
    table_dir.mkdir(parents=True, exist_ok=True)
    html_path.parent.mkdir(parents=True, exist_ok=True)
    md_path.parent.mkdir(parents=True, exist_ok=True)

    if not coloc_brain.empty:
        top_coloc = coloc_brain.sort_values("PP4", ascending=False).head(25).copy()
        top_coloc.insert(0, "rank", range(1, len(top_coloc) + 1))
    else:
        top_coloc = pd.DataFrame(columns=["rank", "locus_id", "tissue", "gene_id", "PP4", "mode"])
    top_coloc.to_csv(table_dir / "top_brain_coloc_hits.tsv", sep="\t", index=False)

    if not mr.empty:
        mr_pres = mr.copy()
        if not mr_rob.empty:
            mr_pres = mr_pres.merge(
                mr_rob[["direction", "cochrans_q_p", "loo_sign_flip_count", "steiger_proxy_supported"]],
                on="direction",
                how="left",
            )
        mr_pres = mr_pres.sort_values(["status", "ivw_p"], na_position="last")
    else:
        mr_pres = pd.DataFrame(columns=["direction", "ivw_beta", "ivw_p", "status", "cochrans_q_p", "loo_sign_flip_count"])
    mr_pres.to_csv(table_dir / "mr_presentation.tsv", sep="\t", index=False)

    sens_rows: list[dict] = []
    if not coloc_sens.empty and {"scenario_id", "n_pp4_ge_0_7", "pp4_max"}.issubset(coloc_sens.columns):
        base = coloc_sens.loc[coloc_sens["scenario_id"].astype(str) == "w500000_p12_1e-05_all"]
        if not base.empty:
            b = base.iloc[0]
            sens_rows.append(
                {
                    "metric": "baseline_coloc_strong",
                    "value": int(b["n_pp4_ge_0_7"]),
                    "details": f"pp4_max={float(b['pp4_max']):.4f}",
                }
            )
        excl = coloc_sens.loc[coloc_sens["mhc_mode"].astype(str) == "exclude_mhc"]
        if not excl.empty:
            e = excl.iloc[0]
            sens_rows.append(
                {
                    "metric": "exclude_mhc_coloc_strong",
                    "value": int(e["n_pp4_ge_0_7"]),
                    "details": f"pp4_max={float(e['pp4_max']):.4f}",
                }
            )
        sens_rows.append(
            {
                "metric": "window_prior_range_coloc_strong",
                "value": int(coloc_sens["n_pp4_ge_0_7"].max() - coloc_sens["n_pp4_ge_0_7"].min()),
                "details": f"min={int(coloc_sens['n_pp4_ge_0_7'].min())}; max={int(coloc_sens['n_pp4_ge_0_7'].max())}",
            }
        )
    if not harm_sens.empty and "n_added_palindromic_window_rows" in harm_sens.columns:
        h = harm_sens.iloc[0]
        sens_rows.append(
            {
                "metric": "added_palindromic_rows",
                "value": int(h["n_added_palindromic_window_rows"]),
                "details": f"top_hit_changed={bool(h.get('top_hit_changed', False))}",
            }
        )
    sensitivity_table = pd.DataFrame(sens_rows)
    sensitivity_table.to_csv(table_dir / "sensitivity_key_metrics.tsv", sep="\t", index=False)

    parts = []
    for title, frame in [
        ("Top brain coloc hits", top_coloc.head(15)),
        ("MR presentation table", mr_pres.head(12)),
        ("Sensitivity key metrics", sensitivity_table),
    ]:
        parts.append(f"<h2>{title}</h2>")
        if frame.empty:
            parts.append("<p>No data</p>")
        else:
            styled = (
                frame.style.hide(axis="index")
                .set_table_styles(
                    [
                        {"selector": "th", "props": "background-color:#111827;color:white;padding:6px;"},
                        {"selector": "td", "props": "padding:6px;border:1px solid #e5e7eb;"},
                    ]
                )
                .to_html()
            )
            parts.append(styled)
    html = "<html><head><meta charset='utf-8'><title>Presentation Tables</title></head><body>" + "".join(parts) + "</body></html>"
    html_path.write_text(html, encoding="utf-8")

    md_lines = ["# Presentation Tables", "", "## Top brain coloc hits", ""]
    md_lines.append(top_coloc.head(15).to_markdown(index=False) if not top_coloc.empty else "No data")
    md_lines.extend(["", "## MR presentation table", ""])
    md_lines.append(mr_pres.head(12).to_markdown(index=False) if not mr_pres.empty else "No data")
    md_lines.extend(["", "## Sensitivity key metrics", ""])
    md_lines.append(sensitivity_table.to_markdown(index=False) if not sensitivity_table.empty else "No data")
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shared", required=True)
    parser.add_argument("--het", required=True)
    parser.add_argument("--coloc-brain", required=True)
    parser.add_argument("--rg", required=True)
    parser.add_argument("--mr", required=True)
    parser.add_argument("--coloc-sensitivity", required=True)
    parser.add_argument("--mr-robustness", required=True)
    parser.add_argument("--harmonization-sensitivity", required=True)
    parser.add_argument("--fig-dir", required=True)
    parser.add_argument("--table-dir", required=True)
    parser.add_argument("--out-html", required=True)
    parser.add_argument("--out-md", required=True)
    args = parser.parse_args()

    style()
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    shared = read_tsv(args.shared)
    het = read_tsv(args.het)
    coloc_brain = read_tsv(args.coloc_brain)
    rg = read_tsv(args.rg)
    mr = read_tsv(args.mr)
    coloc_sens = read_tsv(args.coloc_sensitivity)
    mr_rob = read_tsv(args.mr_robustness)
    harm_sens = read_tsv(args.harmonization_sensitivity)

    fig_upset(shared, fig_dir / "01_colon_overlap_upset.png")
    fig_venn(shared, fig_dir / "02_colon_overlap_venn.png")
    fig_beta_hex(shared, fig_dir / "03_beta_concordance_hexbin.png")
    fig_heterogeneity(het, fig_dir / "04_heterogeneity_volcano_annotated.png")
    fig_coloc_heatmap(coloc_brain, fig_dir / "05_brain_coloc_heatmap.png")
    fig_coloc_locus_bubble(coloc_brain, fig_dir / "06_brain_coloc_locus_bubble.png")
    fig_coloc_sensitivity(coloc_sens, fig_dir / "07_coloc_sensitivity_heatmap.png")
    fig_mhc_sensitivity(coloc_sens, fig_dir / "08_mhc_sensitivity_bar.png")
    fig_rg_panel(rg, fig_dir / "09_genetic_correlation_lollipop.png")
    fig_mr_robustness(mr, mr_rob, fig_dir / "10_mr_forest_robustness.png")

    build_tables(
        coloc_brain=coloc_brain,
        mr=mr,
        mr_rob=mr_rob,
        coloc_sens=coloc_sens,
        harm_sens=harm_sens,
        table_dir=Path(args.table_dir),
        html_path=Path(args.out_html),
        md_path=Path(args.out_md),
    )

    summary = {
        "figures_generated": 10,
        "figure_dir": str(fig_dir),
        "tables_dir": str(Path(args.table_dir)),
    }
    (Path(args.table_dir) / "presentation_manifest.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
