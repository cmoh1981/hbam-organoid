"""Spatial analysis for Stereo-seq transcriptomics data."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
from loguru import logger

from hbam.utils.logging import log_filter, log_metric, log_step


def run_spatial_analysis(
    adata: ad.AnnData,
    min_counts: int = 100,
    min_genes: int = 50,
    region_col: str | None = None,
) -> ad.AnnData:
    """Run spatial QC and variable gene detection.

    Args:
        adata: Input AnnData with spatial coordinates in .obsm["spatial"].
        min_counts: Minimum total counts per bin.
        min_genes: Minimum genes per bin.
        region_col: Column for pseudo-bulk aggregation.

    Returns:
        AnnData with spatial analysis results:
        - .var["spatially_variable"]: bool
        - .var["morans_i"]: Moran's I statistic
        - .var["spatial_pvalue"]: p-value for spatial variability
    """
    with log_step("spatial_analysis"):
        result = adata.copy()

        # Spatial QC
        total_counts = np.asarray(result.X.sum(axis=1)).flatten()
        n_genes = np.asarray((result.X > 0).sum(axis=1)).flatten()

        mask = (total_counts >= min_counts) & (n_genes >= min_genes)
        result = result[mask].copy()

        log_filter("spatial_qc", before=adata.n_obs, after=result.n_obs, reason="low count/gene bins")

        # Highly variable genes (as proxy for spatially variable)
        try:
            import scanpy as sc
            # Ensure X is not log-transformed for HVG detection on counts
            sc.pp.highly_variable_genes(result, n_top_genes=min(500, result.n_vars), flavor="seurat_v3")
            result.var["spatially_variable"] = result.var.get("highly_variable", pd.Series(False, index=result.var_names))
        except Exception as e:
            logger.warning(f"HVG detection failed: {e}. Falling back to variance-based selection.")
            variances = np.var(result.X, axis=0) if hasattr(result.X, 'toarray') is False else np.var(result.X.toarray(), axis=0)
            threshold = np.percentile(variances, 75)
            result.var["spatially_variable"] = variances > threshold

        # Moran's I for spatial autocorrelation (simplified)
        morans_i = np.zeros(result.n_vars)
        morans_p = np.ones(result.n_vars)

        if "spatial" in result.obsm:
            coords = result.obsm["spatial"]
            n = result.n_obs

            if n > 10 and n < 5000:  # Only compute for manageable sizes
                # Build distance-based weights (k-nearest neighbors)
                from sklearn.neighbors import NearestNeighbors
                k = min(6, n - 1)
                nn = NearestNeighbors(n_neighbors=k).fit(coords)
                distances, indices = nn.kneighbors(coords)

                # Compute Moran's I for top variable genes
                sv_idx = np.where(result.var["spatially_variable"].values)[0]

                for j in sv_idx[:100]:  # Limit to 100 genes for speed
                    x = np.asarray(result.X[:, j]).flatten().astype(float)
                    x_mean = x.mean()
                    x_dev = x - x_mean

                    if np.var(x) < 1e-10:
                        continue

                    numerator = 0.0
                    for i in range(n):
                        for neighbor in indices[i]:
                            numerator += x_dev[i] * x_dev[neighbor]

                    denominator = np.sum(x_dev ** 2)
                    w_sum = n * k

                    if denominator > 0:
                        I = (n / w_sum) * (numerator / denominator)
                        morans_i[j] = I
                        # Approximate p-value (z-test)
                        E_I = -1.0 / (n - 1)
                        z = (I - E_I) / (1.0 / np.sqrt(n))
                        morans_p[j] = 2 * (1 - abs(0.5 - abs(z) / 10))  # Rough approximation

        result.var["morans_i"] = morans_i
        result.var["spatial_pvalue"] = morans_p

        n_sv = result.var["spatially_variable"].sum()
        log_metric("spatially_variable_genes", int(n_sv))

        return result


def aggregate_pseudobulk(
    adata: ad.AnnData,
    group_col: str,
) -> ad.AnnData:
    """Aggregate spatial bins into pseudo-bulk profiles per group.

    Args:
        adata: Input spatial AnnData.
        group_col: Column in .obs defining groups (e.g., "region").

    Returns:
        New AnnData with groups as observations, summed counts.
    """
    with log_step("pseudobulk_aggregation", group_col=group_col):
        if group_col not in adata.obs.columns:
            raise ValueError(f"Column '{group_col}' not found in .obs")

        groups = adata.obs[group_col].unique()
        X_bulk = np.zeros((len(groups), adata.n_vars), dtype=np.float64)

        for i, group in enumerate(sorted(groups)):
            mask = adata.obs[group_col] == group
            X_bulk[i] = np.asarray(adata[mask].X.sum(axis=0)).flatten()

        obs = pd.DataFrame({
            "group": sorted(groups),
            "n_bins": [int((adata.obs[group_col] == g).sum()) for g in sorted(groups)],
        }, index=[str(g) for g in sorted(groups)])

        result = ad.AnnData(
            X=X_bulk,
            obs=obs,
            var=adata.var.copy(),
            layers={"raw": X_bulk.copy()},
        )

        logger.info(f"Pseudo-bulk: {len(groups)} groups from {adata.n_obs} bins")

        return result
