"""Missing value imputation with MNAR/MCAR awareness."""

from __future__ import annotations

import anndata as ad
import numpy as np
from loguru import logger

from hbam.utils.logging import log_metric, log_step


def detect_missingness_type(adata: ad.AnnData) -> str:
    """Detect whether missing values are MNAR (left-censored) or MCAR.

    Heuristic: if missing values are enriched in low-intensity features,
    the missingness is likely MNAR (below detection limit).

    Args:
        adata: Input AnnData.

    Returns:
        "MNAR", "MCAR", or "MIXED".
    """
    X = adata.X
    if not np.issubdtype(X.dtype, np.floating):
        return "MCAR"

    is_nan = np.isnan(X)
    n_missing = is_nan.sum()

    if n_missing == 0:
        return "MCAR"  # No missing values

    # Compare mean intensity of features with high vs low missingness
    feature_means = np.nanmean(X, axis=0)
    feature_missingness = is_nan.mean(axis=0)

    # Correlation between feature intensity and missingness
    valid = ~np.isnan(feature_means) & ~np.isnan(feature_missingness)
    if valid.sum() < 10:
        return "MCAR"

    from scipy.stats import spearmanr
    r, p = spearmanr(feature_means[valid], feature_missingness[valid])

    log_metric("missingness_intensity_correlation", f"r={r:.3f}, p={p:.4f}")

    if r < -0.3 and p < 0.05:
        return "MNAR"  # Missing values enriched in low-intensity features
    elif abs(r) < 0.1:
        return "MCAR"
    else:
        return "MIXED"


def impute_mnar(
    adata: ad.AnnData,
    method: str = "minprob",
    q: float = 0.01,
    seed: int = 42,
) -> ad.AnnData:
    """Impute MNAR missing values using left-censored methods.

    MinProb: draw from N(q-th quantile of observed, sd * 0.3).

    Args:
        adata: Input AnnData.
        method: "minprob" or "min".
        q: Quantile for MinProb distribution center.
        seed: Random seed.

    Returns:
        Imputed AnnData with .layers["imputed_mask"].
    """
    with log_step("impute_mnar", method=method, q=q):
        rng = np.random.default_rng(seed)
        result = adata.copy()
        X = result.X.copy()

        is_nan = np.isnan(X)
        n_missing = is_nan.sum()

        if n_missing == 0:
            result.layers["imputed_mask"] = is_nan
            return result

        if method == "minprob":
            for j in range(X.shape[1]):
                col = X[:, j]
                col_nan = np.isnan(col)
                if not col_nan.any():
                    continue

                observed = col[~col_nan]
                if len(observed) < 2:
                    X[col_nan, j] = 0.0
                    continue

                loc = np.quantile(observed, q)
                scale = np.std(observed) * 0.3

                n_fill = col_nan.sum()
                imputed = rng.normal(loc=loc, scale=max(scale, 1e-6), size=n_fill)
                X[col_nan, j] = imputed
        elif method == "min":
            for j in range(X.shape[1]):
                col = X[:, j]
                col_nan = np.isnan(col)
                if col_nan.any():
                    observed = col[~col_nan]
                    X[col_nan, j] = np.min(observed) if len(observed) > 0 else 0.0

        result.X = X
        result.layers["imputed_mask"] = is_nan.astype(np.float32)

        log_metric("imputed_values", int(n_missing))

        return result


def impute_mcar(
    adata: ad.AnnData,
    method: str = "knn",
    k: int = 5,
) -> ad.AnnData:
    """Impute MCAR missing values using KNN or simple methods.

    Args:
        adata: Input AnnData.
        method: "knn", "mean", or "median".
        k: Number of neighbors for KNN.

    Returns:
        Imputed AnnData with .layers["imputed_mask"].
    """
    with log_step("impute_mcar", method=method):
        result = adata.copy()
        X = result.X.copy()

        is_nan = np.isnan(X) if np.issubdtype(X.dtype, np.floating) else np.zeros_like(X, dtype=bool)
        n_missing = is_nan.sum()

        if n_missing == 0:
            result.layers["imputed_mask"] = is_nan.astype(np.float32)
            return result

        if method == "knn":
            from sklearn.impute import KNNImputer
            imputer = KNNImputer(n_neighbors=min(k, X.shape[0] - 1))
            X = imputer.fit_transform(X)
        elif method == "mean":
            col_means = np.nanmean(X, axis=0)
            for j in range(X.shape[1]):
                X[np.isnan(X[:, j]), j] = col_means[j]
        elif method == "median":
            col_medians = np.nanmedian(X, axis=0)
            for j in range(X.shape[1]):
                X[np.isnan(X[:, j]), j] = col_medians[j]

        result.X = X
        result.layers["imputed_mask"] = is_nan.astype(np.float32)

        log_metric("imputed_values", int(n_missing))

        return result


def impute_auto(
    adata: ad.AnnData,
    seed: int = 42,
    mnar_method: str = "minprob",
    mnar_quantile: float = 0.01,
    mcar_method: str = "knn",
    mcar_k: int = 5,
) -> ad.AnnData:
    """Auto-detect missingness type and apply appropriate imputation.

    Args:
        adata: Input AnnData.
        seed: Random seed for MNAR imputation.
        mnar_method: Method for MNAR imputation.
        mnar_quantile: Quantile for MinProb.
        mcar_method: Method for MCAR imputation.
        mcar_k: K for KNN imputation.

    Returns:
        Imputed AnnData with .layers["imputed_mask"].
    """
    with log_step("impute_auto"):
        miss_type = detect_missingness_type(adata)
        logger.info(f"Detected missingness type: {miss_type}")

        if miss_type == "MNAR":
            return impute_mnar(adata, method=mnar_method, q=mnar_quantile, seed=seed)
        elif miss_type == "MCAR":
            return impute_mcar(adata, method=mcar_method, k=mcar_k)
        else:  # MIXED
            logger.info("Mixed missingness: applying MNAR first, then MCAR for remaining")
            result = impute_mnar(adata, method=mnar_method, q=mnar_quantile, seed=seed)
            # Check for remaining NaN
            if np.isnan(result.X).any():
                result = impute_mcar(result, method=mcar_method, k=mcar_k)
            return result
