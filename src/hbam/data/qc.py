"""Quality control filtering for multi-omics data."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
from loguru import logger

from hbam.utils.logging import log_filter, log_metric, log_step


def filter_missingness(
    adata: ad.AnnData,
    threshold: float = 0.7,
    axis: str = "var",
) -> ad.AnnData:
    """Remove features/samples with too many missing values.

    Args:
        adata: Input AnnData.
        threshold: Minimum fraction of non-missing values required.
        axis: "var" to filter features, "obs" to filter samples.

    Returns:
        Filtered AnnData (copy).
    """
    with log_step("filter_missingness", threshold=threshold, axis=axis):
        X = adata.X
        is_nan = np.isnan(X) if np.issubdtype(X.dtype, np.floating) else np.zeros_like(X, dtype=bool)

        if axis == "var":
            presence = 1.0 - is_nan.mean(axis=0)
            mask = presence >= threshold
            result = adata[:, mask].copy()
            log_filter(
                "missingness_var",
                before=adata.n_vars,
                after=result.n_vars,
                reason=f"features present in <{threshold:.0%} of samples",
            )
        elif axis == "obs":
            presence = 1.0 - is_nan.mean(axis=1)
            mask = presence >= threshold
            result = adata[mask].copy()
            log_filter(
                "missingness_obs",
                before=adata.n_obs,
                after=result.n_obs,
                reason=f"samples with <{threshold:.0%} feature coverage",
            )
        else:
            raise ValueError(f"axis must be 'var' or 'obs', got '{axis}'")

        return result


def filter_variance(
    adata: ad.AnnData,
    percentile: float = 10.0,
) -> ad.AnnData:
    """Remove features with near-zero variance.

    Args:
        adata: Input AnnData.
        percentile: Remove features below this variance percentile.

    Returns:
        Filtered AnnData (copy).
    """
    with log_step("filter_variance", percentile=percentile):
        X = adata.X
        variances = np.nanvar(X, axis=0)
        cutoff = np.nanpercentile(variances, percentile)
        mask = variances > cutoff

        result = adata[:, mask].copy()
        log_filter(
            "variance",
            before=adata.n_vars,
            after=result.n_vars,
            reason=f"variance below {percentile}th percentile (cutoff={cutoff:.4f})",
        )

        return result


def filter_samples(
    adata: ad.AnnData,
    min_features: int = 200,
    max_missingness: float = 0.5,
) -> ad.AnnData:
    """Remove low-quality samples.

    Args:
        adata: Input AnnData.
        min_features: Minimum detected features per sample.
        max_missingness: Maximum fraction of missing values per sample.

    Returns:
        Filtered AnnData (copy).
    """
    with log_step("filter_samples", min_features=min_features):
        X = adata.X
        is_nan = np.isnan(X) if np.issubdtype(X.dtype, np.floating) else np.zeros_like(X, dtype=bool)

        n_detected = (~is_nan & (X != 0)).sum(axis=1)
        missingness = is_nan.mean(axis=1)

        mask = (n_detected >= min_features) & (missingness <= max_missingness)

        result = adata[mask].copy()
        log_filter(
            "sample_qc",
            before=adata.n_obs,
            after=result.n_obs,
            reason=f"min_features={min_features}, max_missingness={max_missingness}",
        )

        return result


def qc_report(adata: ad.AnnData) -> dict:
    """Generate QC summary statistics.

    Args:
        adata: Input AnnData.

    Returns:
        Dictionary with QC metrics.
    """
    X = adata.X
    is_nan = np.isnan(X) if np.issubdtype(X.dtype, np.floating) else np.zeros_like(X, dtype=bool)

    report = {
        "n_samples": adata.n_obs,
        "n_features": adata.n_vars,
        "total_elements": X.size,
        "missing_count": int(is_nan.sum()),
        "missing_fraction": float(is_nan.sum() / X.size) if X.size > 0 else 0.0,
        "per_sample_missingness_mean": float(is_nan.mean(axis=1).mean()),
        "per_feature_missingness_mean": float(is_nan.mean(axis=0).mean()),
        "value_range": [float(np.nanmin(X)), float(np.nanmax(X))] if X.size > 0 else [0, 0],
        "median_per_sample": np.nanmedian(X, axis=1).tolist() if X.size > 0 else [],
    }

    log_metric("qc_samples", report["n_samples"])
    log_metric("qc_features", report["n_features"])
    log_metric("qc_missing_fraction", f"{report['missing_fraction']:.4f}")

    return report
