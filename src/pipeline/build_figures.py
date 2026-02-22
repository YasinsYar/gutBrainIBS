from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml


def save_placeholder(path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.5, "No data", ha="center", va="center", fontsize=12)
    ax.set_title(title)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--shared", required=True)
    parser.add_argument("--het", required=True)
    parser.add_argument("--coloc-colon", required=True)
    parser.add_argument("--coloc-brain", required=True)
    parser.add_argument("--rg", required=True)
    parser.add_argument("--mr", required=True)
    parser.add_argument("--pathway", required=True)
    parser.add_argument("--fig-dir", required=True)
    args = parser.parse_args()

    _ = yaml.safe_load(open(args.config, "r", encoding="utf-8"))
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    shared = pd.read_csv(args.shared, sep="\t") if Path(args.shared).exists() else pd.DataFrame()
    if shared.empty:
        save_placeholder(fig_dir / "colon_shared_specific_barplot.png", "Shared vs tissue-specific")
        save_placeholder(fig_dir / "beta_concordance_scatter.png", "Beta concordance")
    else:
        counts = shared["category"].value_counts().rename_axis("category").reset_index(name="count")
        fig, ax = plt.subplots(figsize=(7, 4))
        sns.barplot(data=counts, x="category", y="count", ax=ax)
        ax.set_title("Colon shared vs tissue-specific eQTL")
        fig.tight_layout()
        fig.savefig(fig_dir / "colon_shared_specific_barplot.png", dpi=150)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6, 6))
        scatter_df = shared.dropna(subset=["beta_sigmoid", "beta_transverse"]).copy()
        if len(scatter_df) > 20000:
            scatter_df = scatter_df.sample(20000, random_state=42)
        sns.scatterplot(data=scatter_df, x="beta_sigmoid", y="beta_transverse", s=8, alpha=0.4, ax=ax)
        lim = np.nanmax(np.abs(shared[["beta_sigmoid", "beta_transverse"]].to_numpy())) if len(shared) else 1
        ax.plot([-lim, lim], [-lim, lim], color="black", linestyle="--", linewidth=1)
        ax.set_title("Beta concordance: sigmoid vs transverse")
        fig.tight_layout()
        fig.savefig(fig_dir / "beta_concordance_scatter.png", dpi=150)
        plt.close(fig)

    het = pd.read_csv(args.het, sep="\t") if Path(args.het).exists() else pd.DataFrame()
    if het.empty:
        save_placeholder(fig_dir / "heterogeneity_volcano.png", "Heterogeneity volcano")
    else:
        plot_df = het.copy()
        plot_df["neglog10p"] = -np.log10(np.clip(plot_df["p_het"], 1e-300, 1.0))
        fig, ax = plt.subplots(figsize=(7, 4))
        sns.scatterplot(data=plot_df.sample(min(len(plot_df), 50000), random_state=42), x="z_het", y="neglog10p", s=8, alpha=0.4, ax=ax)
        ax.set_title("Colon heterogeneity")
        fig.tight_layout()
        fig.savefig(fig_dir / "heterogeneity_volcano.png", dpi=150)
        plt.close(fig)

    coloc_colon = pd.read_csv(args.coloc_colon, sep="\t") if Path(args.coloc_colon).exists() else pd.DataFrame()
    coloc_brain = pd.read_csv(args.coloc_brain, sep="\t") if Path(args.coloc_brain).exists() else pd.DataFrame()
    pp = []
    if not coloc_colon.empty and "PP4" in coloc_colon.columns:
        pp.extend([(x, "colon") for x in coloc_colon["PP4"].dropna().tolist()])
    if not coloc_brain.empty and "PP4" in coloc_brain.columns:
        pp.extend([(x, "brain") for x in coloc_brain["PP4"].dropna().tolist()])
    if not pp:
        save_placeholder(fig_dir / "pp4_distribution.png", "PP4 distribution")
    else:
        pp_df = pd.DataFrame(pp, columns=["PP4", "group"])
        fig, ax = plt.subplots(figsize=(7, 4))
        sns.histplot(data=pp_df, x="PP4", hue="group", bins=30, element="step", stat="density", common_norm=False, ax=ax)
        ax.set_title("Colocalization PP4 distributions")
        fig.tight_layout()
        fig.savefig(fig_dir / "pp4_distribution.png", dpi=150)
        plt.close(fig)

    rg = pd.read_csv(args.rg, sep="\t") if Path(args.rg).exists() else pd.DataFrame()
    if rg.empty:
        save_placeholder(fig_dir / "genetic_correlation_panel.png", "Genetic correlation")
    else:
        fig, ax = plt.subplots(figsize=(7, 4))
        sns.barplot(data=rg, x="trait_id", y="rg_proxy", ax=ax)
        ax.axhline(0, color="black", linewidth=1)
        ax.set_title("IBS vs psychiatric trait correlation")
        fig.tight_layout()
        fig.savefig(fig_dir / "genetic_correlation_panel.png", dpi=150)
        plt.close(fig)

    mr = pd.read_csv(args.mr, sep="\t") if Path(args.mr).exists() else pd.DataFrame()
    if mr.empty:
        save_placeholder(fig_dir / "mr_forest.png", "MR forest")
    else:
        plot_df = mr.dropna(subset=["ivw_beta", "ivw_se"]).copy()
        if plot_df.empty:
            save_placeholder(fig_dir / "mr_forest.png", "MR forest")
        else:
            plot_df = plot_df.reset_index(drop=True)
            fig, ax = plt.subplots(figsize=(8, 4))
            y = np.arange(len(plot_df))
            ax.errorbar(plot_df["ivw_beta"], y, xerr=1.96 * plot_df["ivw_se"], fmt="o", color="tab:blue")
            ax.axvline(0, color="black", linestyle="--", linewidth=1)
            ax.set_yticks(y)
            ax.set_yticklabels(plot_df["direction"])
            ax.set_title("MR IVW estimates")
            fig.tight_layout()
            fig.savefig(fig_dir / "mr_forest.png", dpi=150)
            plt.close(fig)

    pathway = pd.read_csv(args.pathway, sep="\t") if Path(args.pathway).exists() else pd.DataFrame()
    if pathway.empty or "q_value" not in pathway.columns:
        save_placeholder(fig_dir / "pathway_enrichment_dotplot.png", "Pathway enrichment")
    else:
        plot_df = pathway.dropna(subset=["q_value"]).head(15)
        if plot_df.empty:
            save_placeholder(fig_dir / "pathway_enrichment_dotplot.png", "Pathway enrichment")
        else:
            fig, ax = plt.subplots(figsize=(8, 4))
            sns.scatterplot(data=plot_df, x="overlap", y="pathway_id", size=-np.log10(np.clip(plot_df["q_value"], 1e-300, 1.0)), hue="q_value", ax=ax)
            ax.set_title("Pathway enrichment")
            fig.tight_layout()
            fig.savefig(fig_dir / "pathway_enrichment_dotplot.png", dpi=150)
            plt.close(fig)


if __name__ == "__main__":
    main()
