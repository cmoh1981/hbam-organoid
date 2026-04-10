"""Region-specific biomarker identification from spatial HBAM scores."""

from __future__ import annotations

import numpy as np
import pandas as pd
from loguru import logger
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

from hbam.utils.logging import log_metric, log_step


def identify_region_biomarkers(
    adata,
    score_col: str = "hbam_score",
    n_regions: int = 4,
    top_n: int = 20,
) -> pd.DataFrame:
    """Identify genes enriched in high-HBAM vs low-HBAM spatial regions.

    Splits spots into quartiles by HBAM score, then runs DE between
    the top and bottom quartiles to find spatially-resolved aging genes.

    Args:
        adata: Spatial AnnData with HBAM scores in .obs[score_col].
        score_col: Column with HBAM scores.
        n_regions: Number of quartile bins (default 4 = quartiles).
        top_n: Number of top genes to return.

    Returns:
        DataFrame with gene, log2fc, pvalue, fdr, high_mean, low_mean, category.
    """
    with log_step("identify_region_biomarkers", n_regions=n_regions):
        if score_col not in adata.obs.columns:
            logger.warning(f"Score column '{score_col}' not found in .obs")
            return pd.DataFrame(columns=["gene", "log2fc", "pvalue", "fdr"])

        scores = adata.obs[score_col].values

        # Split into quartiles
        thresholds = np.percentile(scores, [100 / n_regions, 100 - 100 / n_regions])
        low_mask = scores <= thresholds[0]  # Bottom quartile (young-like)
        high_mask = scores >= thresholds[1]  # Top quartile (aged-like)

        n_low = low_mask.sum()
        n_high = high_mask.sum()

        logger.info(f"Region split: {n_low} low-HBAM spots, {n_high} high-HBAM spots")

        if n_low < 5 or n_high < 5:
            logger.warning("Too few spots in high/low regions for DE analysis")
            return pd.DataFrame(columns=["gene", "log2fc", "pvalue", "fdr"])

        X = adata.X
        n_genes = adata.n_vars

        pvalues = np.ones(n_genes)
        log2fc = np.zeros(n_genes)
        high_means = np.zeros(n_genes)
        low_means = np.zeros(n_genes)

        for j in range(n_genes):
            col_high = np.asarray(X[high_mask, j]).flatten().astype(float)
            col_low = np.asarray(X[low_mask, j]).flatten().astype(float)

            high_means[j] = np.mean(col_high)
            low_means[j] = np.mean(col_low)

            # Log2 fold change (high / low)
            if low_means[j] > 0:
                log2fc[j] = np.log2((high_means[j] + 1e-10) / (low_means[j] + 1e-10))

            if np.std(col_high) > 0 or np.std(col_low) > 0:
                try:
                    _, p = mannwhitneyu(col_high, col_low, alternative="two-sided")
                    pvalues[j] = p
                except ValueError:
                    pass

        # FDR correction
        _, fdr, _, _ = multipletests(pvalues, method="fdr_bh")

        # Build results
        results = pd.DataFrame({
            "gene": adata.var_names.tolist(),
            "log2fc": log2fc,
            "pvalue": pvalues,
            "fdr": fdr,
            "high_hbam_mean": high_means,
            "low_hbam_mean": low_means,
        })

        # Categorize
        results["direction"] = "neutral"
        results.loc[(results["fdr"] < 0.05) & (results["log2fc"] > 0), "direction"] = "high_in_aged"
        results.loc[(results["fdr"] < 0.05) & (results["log2fc"] < 0), "direction"] = "high_in_young"

        # Sort by significance
        results = results.sort_values("fdr")

        n_sig = (results["fdr"] < 0.05).sum()
        n_aged = (results["direction"] == "high_in_aged").sum()
        n_young = (results["direction"] == "high_in_young").sum()

        log_metric("region_de_significant", int(n_sig))
        log_metric("region_de_high_in_aged", int(n_aged))
        log_metric("region_de_high_in_young", int(n_young))

        return results.head(top_n)


def spatial_hbam_summary(adata, score_col: str = "hbam_score") -> dict:
    """Generate summary statistics of spatial HBAM distribution.

    Args:
        adata: Spatial AnnData with HBAM scores.
        score_col: Score column name.

    Returns:
        Dict with spatial HBAM statistics.
    """
    if score_col not in adata.obs.columns:
        return {"error": f"No {score_col} in .obs"}

    scores = adata.obs[score_col].values

    return {
        "n_spots": len(scores),
        "mean": float(np.mean(scores)),
        "std": float(np.std(scores)),
        "min": float(np.min(scores)),
        "max": float(np.max(scores)),
        "q25": float(np.percentile(scores, 25)),
        "median": float(np.percentile(scores, 50)),
        "q75": float(np.percentile(scores, 75)),
        "n_positive": int((scores > 0).sum()),
        "n_negative": int((scores < 0).sum()),
        "fraction_positive": float((scores > 0).mean()),
        "has_spatial": "spatial" in adata.obsm,
    }
