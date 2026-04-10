"""Publication-ready figure generation for the HBAM pipeline."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from loguru import logger

from hbam.output.style import (
    DOUBLE_COL, FULL_PAGE_HEIGHT, PALETTE, SINGLE_COL,
    get_condition_colors, save_figure, set_publication_style,
)
from hbam.utils.logging import log_step


def fig_pipeline_overview(output_dir: Path, **kwargs: Any) -> list[Path]:
    """Figure 1: Pipeline schematic overview."""
    set_publication_style()
    fig, ax = plt.subplots(figsize=(DOUBLE_COL, 3.5))

    boxes = [
        (0.05, 0.6, "Data Layer\n(Ingestion, QC,\nNormalization)", PALETTE["blue"]),
        (0.22, 0.6, "Modality\nAnalysis\n(3 modules)", PALETTE["orange"]),
        (0.39, 0.6, "Integration\n(Orthologs,\nEmbedding)", PALETTE["green"]),
        (0.56, 0.6, "HBAM Engine\n(Weights,\nScoring)", PALETTE["red"]),
        (0.73, 0.6, "Translation\n(Biomarkers,\nClinical)", PALETTE["purple"]),
        (0.90, 0.6, "Output\n(Figures,\nTables)", PALETTE["cyan"]),
    ]

    for x, y, text, color in boxes:
        rect = plt.Rectangle((x - 0.07, y - 0.15), 0.14, 0.35,
                              facecolor=color, alpha=0.3, edgecolor=color, linewidth=1.5,
                              transform=ax.transAxes)
        ax.add_patch(rect)
        ax.text(x, y, text, ha="center", va="center", fontsize=6,
                fontweight="bold", transform=ax.transAxes)

    # Arrows
    for i in range(5):
        x_start = boxes[i][0] + 0.07
        x_end = boxes[i + 1][0] - 0.07
        ax.annotate("", xy=(x_end, 0.6), xytext=(x_start, 0.6),
                    xycoords="axes fraction", textcoords="axes fraction",
                    arrowprops=dict(arrowstyle="->", color="gray", lw=1.5))

    # Input data labels
    inputs = ["Proteomics\n(MaxQuant/DIA-NN)", "Spatial\n(Stereo-seq)", "Functional\n(Muscle)"]
    for i, inp in enumerate(inputs):
        y = 0.15 + i * 0.12
        ax.text(0.05, y, inp, ha="center", va="center", fontsize=5, style="italic",
                transform=ax.transAxes)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("HBAM Multi-Omics Pipeline", fontsize=12, fontweight="bold", pad=10)

    return save_figure(fig, "fig01_pipeline_overview", output_dir)


def fig_qc_summary(mudata: Any, output_dir: Path, **kwargs: Any) -> list[Path]:
    """Figure 2: QC summary (missingness, distributions)."""
    set_publication_style()
    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE_COL, DOUBLE_COL * 0.7))

    for idx, (name, mod) in enumerate(mudata.mod.items()):
        if idx >= 2:
            break
        ax = axes[0, idx]
        X = mod.X

        # Missingness per sample
        if np.issubdtype(X.dtype, np.floating):
            miss = np.isnan(X).mean(axis=1)
        else:
            miss = np.zeros(X.shape[0])
        ax.bar(range(len(miss)), miss, color=PALETTE["blue"], alpha=0.7)
        ax.set_title(f"{name}: Missingness", fontsize=8)
        ax.set_xlabel("Sample")
        ax.set_ylabel("Fraction missing")

    # Distribution boxplot
    ax = axes[1, 0]
    for idx, (name, mod) in enumerate(mudata.mod.items()):
        if idx >= 3:
            break
        medians = np.nanmedian(mod.X, axis=1)
        color = [PALETTE["blue"], PALETTE["orange"], PALETTE["green"]][idx]
        parts = ax.violinplot([medians], positions=[idx], showmeans=True, showmedians=True)
        for pc in parts["bodies"]:
            pc.set_facecolor(color)
            pc.set_alpha(0.5)
    ax.set_xticks(range(len(mudata.mod)))
    ax.set_xticklabels(list(mudata.mod.keys()), rotation=45, ha="right")
    ax.set_title("Per-sample medians by modality", fontsize=8)
    ax.set_ylabel("Median value")

    # Feature variance
    ax = axes[1, 1]
    for idx, (name, mod) in enumerate(mudata.mod.items()):
        if idx >= 3:
            break
        var = np.nanvar(mod.X, axis=0)
        color = [PALETTE["blue"], PALETTE["orange"], PALETTE["green"]][idx]
        ax.hist(np.log10(var + 1e-10), bins=30, alpha=0.5, label=name, color=color)
    ax.set_title("Feature variance distribution", fontsize=8)
    ax.set_xlabel("log10(variance)")
    ax.legend(fontsize=6)

    fig.tight_layout()
    return save_figure(fig, "fig02_qc_summary", output_dir)


def fig_latent_space(mudata: Any, output_dir: Path, **kwargs: Any) -> list[Path]:
    """Figure 3: Latent space visualization (UMAP/PCA)."""
    set_publication_style()

    factor_key = None
    for key in ["pca_factors", "mofa_factors", "X_mofa"]:
        if key in mudata.obsm:
            factor_key = key
            break

    if factor_key is None:
        logger.warning("No latent factors found for visualization")
        fig, ax = plt.subplots(figsize=(SINGLE_COL, SINGLE_COL))
        ax.text(0.5, 0.5, "No embedding available", ha="center", va="center")
        return save_figure(fig, "fig03_latent_space", output_dir)

    factors = mudata.obsm[factor_key]
    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE_COL, SINGLE_COL))

    # Panel A: by condition
    ax = axes[0]
    first_mod = list(mudata.mod.keys())[0]
    if "condition" in mudata.mod[first_mod].obs.columns:
        conditions = mudata.mod[first_mod].obs["condition"].values[:len(factors)]
        colors = get_condition_colors(list(set(conditions)))
        for cond in sorted(set(conditions)):
            mask = conditions == cond
            ax.scatter(factors[mask, 0], factors[mask, 1],
                       c=colors.get(cond, "gray"), label=cond, s=15, alpha=0.7)
        ax.legend(fontsize=5, markerscale=0.8)
    else:
        ax.scatter(factors[:, 0], factors[:, 1], c=PALETTE["blue"], s=15, alpha=0.7)
    ax.set_title("By condition", fontsize=8)
    ax.set_xlabel("Factor 1")
    ax.set_ylabel("Factor 2")

    # Panel B: all samples
    ax = axes[1]
    ax.scatter(factors[:, 0], factors[:, 1], c=PALETTE["blue"], s=15, alpha=0.7)
    ax.set_title("All samples", fontsize=8)
    ax.set_xlabel("Factor 1")

    # Panel C: colored by HBAM score
    ax = axes[2]
    if "hbam_score" in mudata.mod[first_mod].obs.columns:
        scores = mudata.mod[first_mod].obs["hbam_score"].values[:len(factors)]
        sc = ax.scatter(factors[:, 0], factors[:, 1], c=scores, cmap="RdBu_r", s=15, alpha=0.7)
        plt.colorbar(sc, ax=ax, label="HBAM", shrink=0.8)
    else:
        ax.scatter(factors[:, 0], factors[:, 1], c=PALETTE["blue"], s=15, alpha=0.7)
    ax.set_title("By HBAM score", fontsize=8)
    ax.set_xlabel("Factor 1")

    fig.tight_layout()
    return save_figure(fig, "fig03_latent_space", output_dir)


def fig_weight_heatmap(weights: pd.DataFrame, output_dir: Path, top_n: int = 50, **kwargs: Any) -> list[Path]:
    """Figure 4: Gene weight heatmap."""
    set_publication_style()

    df = weights.copy()
    df["abs_weight"] = df["weight"].abs()
    df = df.sort_values("abs_weight", ascending=False).head(top_n)

    fig, ax = plt.subplots(figsize=(SINGLE_COL, max(4, top_n * 0.12)))

    colors = [PALETTE["maturation"] if c == "maturation" else PALETTE["dysfunction"]
              for c in df["category"]]

    ax.barh(range(len(df)), df["weight"], color=colors, alpha=0.8)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["gene"], fontsize=5)
    ax.set_xlabel("Weight")
    ax.set_title(f"Top {len(df)} Gene Weights", fontsize=10)
    ax.invert_yaxis()

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=PALETTE["maturation"], alpha=0.8, label="Maturation"),
        Patch(facecolor=PALETTE["dysfunction"], alpha=0.8, label="Dysfunction"),
    ]
    ax.legend(handles=legend_elements, fontsize=6, loc="lower right")

    fig.tight_layout()
    return save_figure(fig, "fig04_weight_heatmap", output_dir)


def fig_hbam_distribution(scores: pd.DataFrame, output_dir: Path, group_col: str = "condition", **kwargs: Any) -> list[Path]:
    """Figure 5: HBAM score distribution by condition."""
    set_publication_style()
    fig, ax = plt.subplots(figsize=(SINGLE_COL, SINGLE_COL))

    if group_col in scores.columns:
        conditions = sorted(scores[group_col].unique())
        colors = get_condition_colors(conditions)

        data_groups = [scores[scores[group_col] == c]["hbam_score"].values for c in conditions]
        parts = ax.violinplot(data_groups, showmeans=True, showmedians=True)

        for i, pc in enumerate(parts["bodies"]):
            pc.set_facecolor(colors.get(conditions[i], "gray"))
            pc.set_alpha(0.5)

        # Overlay strip plot
        for i, (cond, data) in enumerate(zip(conditions, data_groups)):
            jitter = np.random.default_rng(42).uniform(-0.1, 0.1, len(data))
            ax.scatter(np.full(len(data), i + 1) + jitter, data,
                       c=colors.get(cond, "gray"), s=10, alpha=0.6, zorder=5)

        ax.set_xticks(range(1, len(conditions) + 1))
        ax.set_xticklabels(conditions)

        # Add p-value annotation
        if len(data_groups) >= 2 and len(data_groups[0]) >= 2 and len(data_groups[1]) >= 2:
            from scipy.stats import mannwhitneyu
            _, p = mannwhitneyu(data_groups[0], data_groups[1], alternative="two-sided")
            y_max = max(max(d) for d in data_groups if len(d) > 0)
            ax.text(1.5, y_max * 1.05, f"p = {p:.2e}", ha="center", fontsize=7)
    else:
        ax.hist(scores["hbam_score"], bins=20, color=PALETTE["blue"], alpha=0.7)

    ax.set_ylabel("HBAM Score")
    ax.set_title("HBAM Score Distribution", fontsize=10)
    fig.tight_layout()
    return save_figure(fig, "fig05_hbam_distribution", output_dir)


def fig_pathway_enrichment(enrichment_results: pd.DataFrame | None, output_dir: Path, top_n: int = 20, **kwargs: Any) -> list[Path]:
    """Figure 6: Pathway enrichment dot plot."""
    set_publication_style()
    fig, ax = plt.subplots(figsize=(DOUBLE_COL, SINGLE_COL * 1.5))

    if enrichment_results is not None and len(enrichment_results) > 0:
        df = enrichment_results.head(top_n).copy()

        term_col = next((c for c in df.columns if "term" in c.lower() or "name" in c.lower()), df.columns[0])
        pval_col = next((c for c in df.columns if "p" in c.lower() and "adj" in c.lower()),
                        next((c for c in df.columns if "fdr" in c.lower()),
                             next((c for c in df.columns if "p" in c.lower()), None)))

        if pval_col:
            df["neg_log_p"] = -np.log10(df[pval_col].clip(1e-20))
            ax.barh(range(len(df)), df["neg_log_p"], color=PALETTE["red"], alpha=0.7)
            ax.set_yticks(range(len(df)))
            ax.set_yticklabels(df[term_col], fontsize=5)
            ax.set_xlabel("-log\u2081\u2080(FDR)")
            ax.axvline(x=-np.log10(0.05), color="gray", linestyle="--", linewidth=0.5, label="FDR=0.05")
        else:
            ax.text(0.5, 0.5, "No enrichment p-values", ha="center", va="center", transform=ax.transAxes)
    else:
        ax.text(0.5, 0.5, "No enrichment results available", ha="center", va="center", transform=ax.transAxes)

    ax.set_title("Pathway Enrichment", fontsize=10)
    ax.invert_yaxis()
    fig.tight_layout()
    return save_figure(fig, "fig06_pathway_enrichment", output_dir)


def fig_volcano(de_results: pd.DataFrame | None, output_dir: Path, fdr_threshold: float = 0.05, fc_threshold: float = 1.0, top_n: int = 10, **kwargs: Any) -> list[Path]:
    """Figure 7: Volcano plot."""
    set_publication_style()
    fig, ax = plt.subplots(figsize=(SINGLE_COL, SINGLE_COL))

    if de_results is not None and len(de_results) > 0:
        df = de_results.copy()

        fc_col = next((c for c in df.columns if "log2fc" in c.lower() or "fold" in c.lower()), None)
        p_col = next((c for c in df.columns if "fdr" in c.lower() or "padj" in c.lower()),
                     next((c for c in df.columns if "pvalue" in c.lower() or "p_value" in c.lower()), None))

        if fc_col and p_col:
            df["neg_log_p"] = -np.log10(df[p_col].clip(1e-300))

            sig_up = (df[p_col] < fdr_threshold) & (df[fc_col] > fc_threshold)
            sig_down = (df[p_col] < fdr_threshold) & (df[fc_col] < -fc_threshold)
            non_sig = ~(sig_up | sig_down)

            ax.scatter(df.loc[non_sig, fc_col], df.loc[non_sig, "neg_log_p"],
                       c=PALETTE["non_significant"], s=5, alpha=0.3)
            ax.scatter(df.loc[sig_up, fc_col], df.loc[sig_up, "neg_log_p"],
                       c=PALETTE["up"], s=8, alpha=0.7, label="Up")
            ax.scatter(df.loc[sig_down, fc_col], df.loc[sig_down, "neg_log_p"],
                       c=PALETTE["down"], s=8, alpha=0.7, label="Down")

            # Label top genes
            sig = df[sig_up | sig_down].nlargest(top_n, "neg_log_p")
            gene_col = next((c for c in df.columns if "gene" in c.lower()), df.index.name or "index")
            for _, row in sig.iterrows():
                gene_name = row.get(gene_col, row.name) if gene_col != "index" else row.name
                ax.annotate(str(gene_name), (row[fc_col], row["neg_log_p"]),
                            fontsize=4, ha="center", va="bottom")

            ax.axhline(-np.log10(fdr_threshold), color="gray", linestyle="--", linewidth=0.5)
            ax.axvline(fc_threshold, color="gray", linestyle="--", linewidth=0.5)
            ax.axvline(-fc_threshold, color="gray", linestyle="--", linewidth=0.5)

            ax.set_xlabel("log\u2082 fold change")
            ax.set_ylabel("-log\u2081\u2080(FDR)")
            ax.legend(fontsize=6)
    else:
        ax.text(0.5, 0.5, "No DE results available", ha="center", va="center", transform=ax.transAxes)

    ax.set_title("Volcano Plot", fontsize=10)
    fig.tight_layout()
    return save_figure(fig, "fig07_volcano", output_dir)


def fig_correlation_matrix(scores: pd.DataFrame, output_dir: Path, **kwargs: Any) -> list[Path]:
    """Figure 8: Correlation matrix of HBAM components."""
    set_publication_style()

    cols = [c for c in ["hbam_score", "dysfunction_score", "maturation_score"] if c in scores.columns]
    if not cols:
        fig, ax = plt.subplots(figsize=(SINGLE_COL, SINGLE_COL))
        ax.text(0.5, 0.5, "No score data", ha="center", va="center", transform=ax.transAxes)
        return save_figure(fig, "fig08_correlation_matrix", output_dir)

    corr = scores[cols].corr()
    fig, ax = plt.subplots(figsize=(SINGLE_COL, SINGLE_COL))
    sns.heatmap(corr, annot=True, fmt=".2f", cmap="RdBu_r", center=0,
                square=True, ax=ax, cbar_kws={"shrink": 0.8},
                xticklabels=[c.replace("_score", "").replace("_", "\n") for c in cols],
                yticklabels=[c.replace("_score", "").replace("_", "\n") for c in cols])
    ax.set_title("Score Correlations", fontsize=10)
    fig.tight_layout()
    return save_figure(fig, "fig08_correlation_matrix", output_dir)


def fig_biomarker_panel(panel: pd.DataFrame, mudata: Any, output_dir: Path, **kwargs: Any) -> list[Path]:
    """Figure 9: Biomarker panel heatmap."""
    set_publication_style()

    first_mod = list(mudata.mod.keys())[0]
    mod = mudata.mod[first_mod]

    panel_genes = [g for g in panel["gene"] if g in mod.var_names]
    if not panel_genes:
        fig, ax = plt.subplots(figsize=(SINGLE_COL, SINGLE_COL))
        ax.text(0.5, 0.5, "No panel genes found", ha="center", va="center", transform=ax.transAxes)
        return save_figure(fig, "fig09_biomarker_panel", output_dir)

    X = mod[:, panel_genes].X.copy()
    if np.isnan(X).any():
        X = np.nan_to_num(X, nan=0.0)

    fig, ax = plt.subplots(figsize=(DOUBLE_COL, max(3, len(panel_genes) * 0.15)))
    sns.heatmap(X.T, cmap="RdBu_r", center=0, ax=ax,
                yticklabels=panel_genes, xticklabels=False,
                cbar_kws={"shrink": 0.5, "label": "Expression"})
    ax.set_title("Biomarker Panel", fontsize=10)
    ax.set_xlabel("Samples")
    fig.tight_layout()
    return save_figure(fig, "fig09_biomarker_panel", output_dir)


def fig_integration_diagnostic(mudata: Any, output_dir: Path, **kwargs: Any) -> list[Path]:
    """Figure 10: Integration diagnostic."""
    set_publication_style()
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE_COL, SINGLE_COL))

    # Variance explained
    ax = axes[0]
    for key in ["pca_variance_explained", "mofa_variance_explained"]:
        if key in mudata.uns:
            var_exp = mudata.uns[key]
            ax.bar(range(len(var_exp)), var_exp, color=PALETTE["blue"], alpha=0.7)
            ax.set_xlabel("Factor/Component")
            ax.set_ylabel("Variance explained")
            ax.set_title("Variance per factor", fontsize=8)
            break
    else:
        ax.text(0.5, 0.5, "No variance data", ha="center", va="center", transform=ax.transAxes)

    # Gene count per modality
    ax = axes[1]
    mod_names = list(mudata.mod.keys())
    n_genes = [mudata.mod[m].n_vars for m in mod_names]
    colors = [PALETTE["blue"], PALETTE["orange"], PALETTE["green"]][:len(mod_names)]
    ax.bar(range(len(mod_names)), n_genes, color=colors, alpha=0.7)
    ax.set_xticks(range(len(mod_names)))
    ax.set_xticklabels(mod_names, rotation=45, ha="right")
    ax.set_ylabel("Number of genes")
    ax.set_title("Genes per modality", fontsize=8)

    fig.tight_layout()
    return save_figure(fig, "fig10_integration_diagnostic", output_dir)
