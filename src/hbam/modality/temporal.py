"""Temporal analysis for liver proteomics data."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
from loguru import logger
from scipy import stats
from statsmodels.stats.multitest import multipletests

from hbam.utils.logging import log_metric, log_step


def run_temporal_analysis(
    adata: ad.AnnData,
    time_col: str = "timepoint",
    condition_col: str | None = None,
) -> ad.AnnData:
    """Run temporal differential expression analysis on proteomics data.

    Performs Kruskal-Wallis test across timepoints and Spearman trend detection.

    Args:
        adata: Input AnnData with temporal metadata in .obs.
        time_col: Column name in .obs containing timepoint info.
        condition_col: Optional condition column for stratified analysis.

    Returns:
        AnnData with DE results in .var columns:
        - temporal_pvalue: raw p-value from Kruskal-Wallis
        - temporal_fdr: BH-adjusted p-value
        - temporal_fold_change: max fold change across timepoints
        - temporal_trend: trend direction ("up", "down", "flat")
        - temporal_spearman_r: Spearman correlation with time
    """
    with log_step("temporal_analysis", time_col=time_col):
        result = adata.copy()

        if time_col not in result.obs.columns:
            logger.warning(f"Time column '{time_col}' not found in .obs. Using index as proxy.")
            result.obs[time_col] = range(result.n_obs)

        timepoints = result.obs[time_col].values
        unique_times = np.unique(timepoints)
        n_timepoints = len(unique_times)

        logger.info(f"Temporal analysis: {result.n_vars} features, {n_timepoints} timepoints")

        pvalues = np.ones(result.n_vars)
        fold_changes = np.zeros(result.n_vars)
        spearman_r = np.zeros(result.n_vars)
        trends = np.full(result.n_vars, "flat", dtype=object)

        X = result.X
        time_numeric = pd.to_numeric(timepoints, errors="coerce")
        if np.any(np.isnan(time_numeric)):
            # Map to ordinal values
            time_map = {t: i for i, t in enumerate(sorted(unique_times))}
            time_numeric = np.array([time_map[t] for t in timepoints], dtype=float)

        for j in range(result.n_vars):
            col = X[:, j]
            valid = ~np.isnan(col) if np.issubdtype(col.dtype, np.floating) else np.ones(len(col), dtype=bool)

            if valid.sum() < 3:
                continue

            col_valid = col[valid]
            time_valid = time_numeric[valid]

            # Kruskal-Wallis across timepoints
            if n_timepoints >= 2:
                groups = []
                for t in unique_times:
                    t_idx = timepoints[valid] == t
                    if t_idx.sum() > 0:
                        groups.append(col_valid[t_idx])

                if len(groups) >= 2 and all(len(g) > 0 for g in groups):
                    try:
                        stat, p = stats.kruskal(*groups)
                        pvalues[j] = p
                    except ValueError:
                        pass

            # Fold change (max group mean / min group mean)
            group_means = [np.nanmean(g) for g in groups if len(g) > 0]
            if group_means and min(group_means) > 0:
                fold_changes[j] = max(group_means) / min(group_means)

            # Trend detection via Spearman correlation
            if len(col_valid) >= 3:
                r, p_trend = stats.spearmanr(time_valid, col_valid)
                spearman_r[j] = r if np.isfinite(r) else 0.0

                if abs(r) > 0.3 and p_trend < 0.05:
                    trends[j] = "up" if r > 0 else "down"

        # FDR correction
        reject, fdr, _, _ = multipletests(pvalues, method="fdr_bh")

        result.var["temporal_pvalue"] = pvalues
        result.var["temporal_fdr"] = fdr
        result.var["temporal_fold_change"] = fold_changes
        result.var["temporal_trend"] = trends
        result.var["temporal_spearman_r"] = spearman_r

        n_sig = (fdr < 0.05).sum()
        n_up = (trends == "up").sum()
        n_down = (trends == "down").sum()

        log_metric("temporal_significant_features", int(n_sig), fdr_threshold=0.05)
        log_metric("temporal_trends", f"up={n_up}, down={n_down}, flat={result.n_vars - n_up - n_down}")

        return result
