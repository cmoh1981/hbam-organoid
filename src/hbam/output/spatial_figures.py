"""Spatial HBAM visualization figures."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from loguru import logger

from hbam.output.style import (
    DOUBLE_COL, PALETTE, SINGLE_COL,
    save_figure, set_publication_style,
)
from hbam.utils.logging import log_step


def fig_spatial_hbam_map(
    adata: Any,
    output_dir: Path,
    score_col: str = "hbam_score",
    **kwargs: Any,
) -> list[Path]:
    """Spatial scatter plot of HBAM scores on tissue coordinates.

    4 panels: (A) HBAM score, (B) Dysfunction, (C) Maturation, (D) Top gene.

    Args:
        adata: Spatial AnnData with scores in .obs and coords in .obsm["spatial"].
        output_dir: Output directory for figures.
        score_col: Column name for HBAM score in .obs.

    Returns:
        List of saved file paths.
    """
    set_publication_style()

    if "spatial" not in adata.obsm:
        logger.warning("No spatial coordinates found. Skipping spatial HBAM map.")
        fig, ax = plt.subplots(figsize=(SINGLE_COL, SINGLE_COL))
        ax.text(0.5, 0.5, "No spatial data", ha="center", va="center", transform=ax.transAxes)
        return save_figure(fig, "fig_spatial_hbam_map", output_dir)

    coords = adata.obsm["spatial"]
    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE_COL, DOUBLE_COL))

    score_cols = [
        ("hbam_score", "HBAM Score\n(aging burden)", "RdBu_r"),
        ("dysfunction_score", "Dysfunction Score", "Reds"),
        ("maturation_score", "Maturation Score", "Blues"),
    ]

    for idx, (col, title, cmap) in enumerate(score_cols):
        ax = axes.flat[idx]
        if col in adata.obs.columns:
            scores = adata.obs[col].values
            # Clip outliers for better visualization
            vmin, vmax = np.percentile(scores, [2, 98])
            sc = ax.scatter(
                coords[:, 0], coords[:, 1],
                c=scores, cmap=cmap, s=1, alpha=0.8,
                vmin=vmin, vmax=vmax, rasterized=True,
            )
            plt.colorbar(sc, ax=ax, shrink=0.7, label=col.replace("_", " ").title())
        else:
            ax.text(0.5, 0.5, f"No {col}", ha="center", va="center", transform=ax.transAxes)

        ax.set_title(title, fontsize=8)
        ax.set_aspect("equal")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")

    # Panel D: spot density or n_genes
    ax = axes[1, 1]
    if "n_genes_detected" in adata.obs.columns:
        vals = adata.obs["n_genes_detected"].values
    else:
        vals = np.asarray((adata.X > 0).sum(axis=1)).flatten()

    sc = ax.scatter(
        coords[:, 0], coords[:, 1],
        c=vals, cmap="viridis", s=1, alpha=0.8, rasterized=True,
    )
    plt.colorbar(sc, ax=ax, shrink=0.7, label="Genes detected")
    ax.set_title("Gene Coverage", fontsize=8)
    ax.set_aspect("equal")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    fig.suptitle("Spatial HBAM Analysis", fontsize=12, fontweight="bold", y=1.02)
    fig.tight_layout()
    return save_figure(fig, "fig_spatial_hbam_map", output_dir)


def fig_spatial_gene_overlay(
    adata: Any,
    genes: list[str],
    output_dir: Path,
    n_cols: int = 3,
    **kwargs: Any,
) -> list[Path]:
    """Overlay top aging gene expression on spatial coordinates.

    Args:
        adata: Spatial AnnData.
        genes: List of genes to visualize.
        output_dir: Output directory.
        n_cols: Number of columns in grid.

    Returns:
        List of saved file paths.
    """
    set_publication_style()

    if "spatial" not in adata.obsm:
        logger.warning("No spatial coordinates. Skipping gene overlay.")
        fig, ax = plt.subplots(figsize=(SINGLE_COL, SINGLE_COL))
        ax.text(0.5, 0.5, "No spatial data", ha="center", va="center", transform=ax.transAxes)
        return save_figure(fig, "fig_spatial_gene_overlay", output_dir)

    # Filter to genes that exist in the data
    available = [g for g in genes if g in adata.var_names]
    if not available:
        # Try case-insensitive
        var_upper = {v.upper(): v for v in adata.var_names}
        available = [var_upper[g.upper()] for g in genes if g.upper() in var_upper]

    n_genes = min(len(available), 6)
    if n_genes == 0:
        logger.warning("No matching genes found for spatial overlay.")
        fig, ax = plt.subplots(figsize=(SINGLE_COL, SINGLE_COL))
        ax.text(0.5, 0.5, "No matching genes", ha="center", va="center", transform=ax.transAxes)
        return save_figure(fig, "fig_spatial_gene_overlay", output_dir)

    n_rows = (n_genes + n_cols - 1) // n_cols
    coords = adata.obsm["spatial"]

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(DOUBLE_COL, SINGLE_COL * n_rows))
    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)

    for idx in range(n_rows * n_cols):
        row, col = idx // n_cols, idx % n_cols
        ax = axes[row, col]

        if idx < n_genes:
            gene = available[idx]
            gene_idx = list(adata.var_names).index(gene)
            expr = np.asarray(adata.X[:, gene_idx]).flatten().astype(float)

            vmin, vmax = np.percentile(expr[expr > 0], [5, 95]) if (expr > 0).any() else (0, 1)

            sc = ax.scatter(
                coords[:, 0], coords[:, 1],
                c=expr, cmap="magma", s=1, alpha=0.8,
                vmin=max(0, vmin), vmax=vmax, rasterized=True,
            )
            plt.colorbar(sc, ax=ax, shrink=0.7)
            ax.set_title(gene, fontsize=8, fontweight="bold")
        else:
            ax.axis("off")

        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle("Top Aging Genes — Spatial Expression", fontsize=10, fontweight="bold")
    fig.tight_layout()
    return save_figure(fig, "fig_spatial_gene_overlay", output_dir)
