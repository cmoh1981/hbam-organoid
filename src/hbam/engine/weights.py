"""Gene weight computation from latent factor loadings."""

from __future__ import annotations

import mudata as md
import numpy as np
import pandas as pd
from loguru import logger

from hbam.utils.logging import log_metric, log_step


def compute_weights_mofa(
    mudata: md.MuData,
    gene_sets: dict[str, set[str]],
    factors: list[int] | None = None,
) -> pd.DataFrame:
    """Compute gene weights from MOFA+ factor loadings.

    Args:
        mudata: MuData with MOFA+ results.
        gene_sets: Dict with "maturation" and "dysfunction" gene sets.
        factors: Which factors to use. None = auto-select.

    Returns:
        DataFrame with columns: gene, weight, factor, category, modality.
    """
    with log_step("compute_weights_mofa"):
        rows = []

        for mod_name, mod in mudata.mod.items():
            loadings_key = None
            for key in ["mofa_weights", "LFs", "pca_loadings"]:
                if key in mod.varm:
                    loadings_key = key
                    break

            if loadings_key is None:
                logger.warning(f"No loadings found for modality '{mod_name}', skipping")
                continue

            loadings = mod.varm[loadings_key]

            if factors is None:
                use_factors = list(range(min(5, loadings.shape[1])))
            else:
                use_factors = factors

            for j, gene in enumerate(mod.var_names):
                gene_upper = gene.upper()

                if gene_upper in gene_sets.get("maturation", set()):
                    category = "maturation"
                elif gene_upper in gene_sets.get("dysfunction", set()):
                    category = "dysfunction"
                else:
                    continue

                weight = np.mean(np.abs(loadings[j, use_factors]))

                rows.append({
                    "gene": gene,
                    "weight": float(weight),
                    "category": category,
                    "modality": mod_name,
                    "factor_indices": str(use_factors),
                })

        df = pd.DataFrame(rows)

        if len(df) == 0:
            logger.warning("No genes matched gene sets. Returning empty weights.")
            return pd.DataFrame(columns=["gene", "weight", "category", "modality", "factor_indices"])

        log_metric("weighted_genes", len(df))
        log_metric("weighted_maturation", (df["category"] == "maturation").sum())
        log_metric("weighted_dysfunction", (df["category"] == "dysfunction").sum())

        return df


def compute_weights_pca(
    mudata: md.MuData,
    gene_sets: dict[str, set[str]],
    components: list[int] | None = None,
) -> pd.DataFrame:
    """Compute gene weights from PCA loadings.

    Args:
        mudata: MuData with PCA results.
        gene_sets: Dict with "maturation" and "dysfunction" gene sets.
        components: Which components to use.

    Returns:
        DataFrame with columns: gene, weight, category, modality.
    """
    with log_step("compute_weights_pca"):
        rows = []

        for mod_name, mod in mudata.mod.items():
            if "pca_loadings" not in mod.varm:
                logger.warning(f"No PCA loadings for '{mod_name}', skipping")
                continue

            loadings = mod.varm["pca_loadings"]

            if components is None:
                use_comps = list(range(min(5, loadings.shape[1])))
            else:
                use_comps = components

            for j, gene in enumerate(mod.var_names):
                gene_upper = gene.upper()

                if gene_upper in gene_sets.get("maturation", set()):
                    category = "maturation"
                elif gene_upper in gene_sets.get("dysfunction", set()):
                    category = "dysfunction"
                else:
                    continue

                weight = np.mean(np.abs(loadings[j, use_comps]))

                rows.append({
                    "gene": gene,
                    "weight": float(weight),
                    "category": category,
                    "modality": mod_name,
                })

        return pd.DataFrame(rows) if rows else pd.DataFrame(columns=["gene", "weight", "category", "modality"])


def normalize_weights(weights: pd.DataFrame) -> pd.DataFrame:
    """Normalize weights to sum to 1 within each category.

    Args:
        weights: DataFrame with gene, weight, category columns.

    Returns:
        DataFrame with normalized weights.
    """
    result = weights.copy()

    for cat in result["category"].unique():
        mask = result["category"] == cat
        total = result.loc[mask, "weight"].sum()
        if total > 0:
            result.loc[mask, "weight"] = result.loc[mask, "weight"] / total

    return result
