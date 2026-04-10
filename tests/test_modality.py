"""Tests for modality analysis modules."""

import numpy as np
import pytest
import sys
sys.path.insert(0, "src")


def test_temporal_analysis(sample_proteomics):
    from hbam.modality.temporal import run_temporal_analysis

    result = run_temporal_analysis(sample_proteomics, time_col="timepoint")

    assert "temporal_pvalue" in result.var.columns
    assert "temporal_fdr" in result.var.columns
    assert "temporal_fold_change" in result.var.columns
    assert "temporal_trend" in result.var.columns
    assert "temporal_spearman_r" in result.var.columns

    # FDR >= raw p-value (BH correction never makes p-values smaller)
    valid = result.var["temporal_pvalue"] < 1.0
    if valid.any():
        assert (
            result.var.loc[valid, "temporal_fdr"]
            >= result.var.loc[valid, "temporal_pvalue"] - 1e-10
        ).all()

    # Trends are valid values
    assert all(t in ("up", "down", "flat") for t in result.var["temporal_trend"])


def test_temporal_no_time_col(sample_proteomics):
    """Test graceful handling when time column is missing."""
    from hbam.modality.temporal import run_temporal_analysis

    # Missing time_col falls back to index-as-proxy, still produces results
    result = run_temporal_analysis(sample_proteomics, time_col="nonexistent")
    assert "temporal_pvalue" in result.var.columns


def test_spatial_analysis(sample_spatial):
    from hbam.modality.spatial import run_spatial_analysis

    result = run_spatial_analysis(sample_spatial, min_counts=1, min_genes=1)

    assert "spatially_variable" in result.var.columns
    assert "morans_i" in result.var.columns
    assert "spatial_pvalue" in result.var.columns

    # Spatial coords preserved
    assert "spatial" in result.obsm

    # No imputation of zeros — zero values should still be present
    assert (result.X == 0).any()


def test_pseudobulk(sample_spatial):
    from hbam.modality.spatial import aggregate_pseudobulk

    result = aggregate_pseudobulk(sample_spatial, group_col="region")

    assert result.n_obs == sample_spatial.obs["region"].nunique()
    assert result.n_vars == sample_spatial.n_vars
    assert "n_bins" in result.obs.columns
    assert "raw" in result.layers


def test_functional_analysis(sample_proteomics):
    from hbam.modality.functional import run_functional_analysis

    result = run_functional_analysis(sample_proteomics, condition_col="condition")

    assert "de_pvalue" in result.var.columns
    assert "de_fdr" in result.var.columns
    assert "log2fc" in result.var.columns
    assert "functional_category" in result.var.columns

    # p-values in valid range
    assert (result.var["de_pvalue"] >= 0).all()
    assert (result.var["de_pvalue"] <= 1).all()
    assert (result.var["de_fdr"] >= 0).all()


def test_functional_no_condition(sample_proteomics):
    """Test graceful handling when condition column is missing."""
    from hbam.modality.functional import run_functional_analysis

    result = run_functional_analysis(sample_proteomics, condition_col="nonexistent")
    assert "de_pvalue" in result.var.columns
    # Should default to 1.0 when condition column is absent
    assert (result.var["de_pvalue"] == 1.0).all()
