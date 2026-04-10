"""Normalization methods for multi-omics data."""

from __future__ import annotations

import anndata as ad
import numpy as np
from loguru import logger

from hbam.utils.logging import log_metric, log_step


def normalize_proteomics(
    adata: ad.AnnData,
    method: str = "median",
    pseudocount: float = 1.0,
) -> ad.AnnData:
    """Normalize proteomics data: median centering then log2 transform.

    Args:
        adata: Input AnnData with raw intensity values.
        method: Normalization method ("median" or "quantile").
        pseudocount: Added before log2 transform.

    Returns:
        Normalized AnnData with .layers["raw"] preserved.
    """
    with log_step("normalize_proteomics", method=method):
        result = adata.copy()
        X = result.X.copy()

        # Store raw if not already stored
        if "raw" not in result.layers:
            result.layers["raw"] = X.copy()

        if method == "median":
            # Median centering in linear space (per sample)
            medians = np.nanmedian(X, axis=1, keepdims=True)
            medians[medians == 0] = 1.0  # Avoid division by zero
            global_median = np.nanmedian(medians)
            X = X / medians * global_median
        elif method == "quantile":
            # Simple quantile normalization
            X = _quantile_normalize(X)
        else:
            logger.warning(f"Unknown normalization method '{method}', skipping normalization")

        # Log2 transform
        X = np.log2(X + pseudocount)

        result.X = X

        # Report metrics
        sample_medians = np.nanmedian(result.X, axis=1)
        log_metric("norm_median_range", f"[{np.min(sample_medians):.2f}, {np.max(sample_medians):.2f}]")

        return result


def normalize_spatial(
    adata: ad.AnnData,
    method: str = "scran",
    target_sum: float = 1e4,
) -> ad.AnnData:
    """Normalize spatial transcriptomics data.

    Args:
        adata: Input AnnData with count data.
        method: "sctransform", "scran", or "total".
        target_sum: Target sum for total-count normalization.

    Returns:
        Normalized AnnData.
    """
    with log_step("normalize_spatial", method=method):
        result = adata.copy()

        if "raw" not in result.layers:
            result.layers["raw"] = result.X.copy()

        if method == "sctransform":
            try:
                import scanpy as sc
                sc.pp.normalize_total(result, target_sum=target_sum)
                sc.pp.log1p(result)
                logger.info("SCTransform fallback: using normalize_total + log1p")
            except ImportError:
                logger.warning("scanpy not available, falling back to total normalization")
                result = _normalize_total(result, target_sum)
        elif method == "scran":
            try:
                import scanpy as sc
                sc.pp.normalize_total(result, target_sum=target_sum)
                sc.pp.log1p(result)
            except ImportError:
                result = _normalize_total(result, target_sum)
        elif method == "total":
            result = _normalize_total(result, target_sum)
        else:
            logger.warning(f"Unknown method '{method}', using total normalization")
            result = _normalize_total(result, target_sum)

        return result


def normalize_generic(
    adata: ad.AnnData,
    method: str = "log2",
    pseudocount: float = 1.0,
) -> ad.AnnData:
    """Generic normalization (log2 transform).

    Args:
        adata: Input AnnData.
        method: "log2" or "log1p".
        pseudocount: Added before log transform.

    Returns:
        Normalized AnnData.
    """
    with log_step("normalize_generic", method=method):
        result = adata.copy()

        if "raw" not in result.layers:
            result.layers["raw"] = result.X.copy()

        if method == "log2":
            result.X = np.log2(result.X + pseudocount)
        elif method == "log1p":
            result.X = np.log1p(result.X)

        return result


def _quantile_normalize(X: np.ndarray) -> np.ndarray:
    """Simple quantile normalization."""
    X_sorted = np.sort(X, axis=0)
    means = np.nanmean(X_sorted, axis=1, keepdims=True)
    ranks = np.argsort(np.argsort(X, axis=0), axis=0)
    result = np.zeros_like(X)
    for j in range(X.shape[1]):
        result[:, j] = means[ranks[:, j], 0]
    return result


def _normalize_total(adata: ad.AnnData, target_sum: float) -> ad.AnnData:
    """Manual total-count normalization + log1p."""
    X = adata.X.astype(np.float64)
    totals = X.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1.0
    X = X / totals * target_sum
    X = np.log1p(X)
    adata.X = X.astype(np.float32)
    return adata
