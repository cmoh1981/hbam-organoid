"""Tests for spatial data integration."""

import numpy as np
import pandas as pd
import pytest
import sys
sys.path.insert(0, "src")


@pytest.fixture
def spatial_adata_with_scores():
    """Create spatial AnnData with HBAM scores for testing."""
    import anndata as ad

    rng = np.random.default_rng(42)
    n_spots = 200
    n_genes = 100

    X = rng.poisson(5, size=(n_spots, n_genes)).astype(np.float32)
    coords = np.column_stack([
        rng.uniform(0, 100, n_spots),
        rng.uniform(0, 100, n_spots),
    ]).astype(np.float32)

    gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]
    # Add some real gene names
    real = ["ITIH2", "SELENBP1", "C3", "ACLY", "FASN", "TXNRD1", "ETFA", "NDUFB4", "CD63", "MFAP4"]
    for i, g in enumerate(real[:min(len(real), n_genes)]):
        gene_names[i] = g

    obs = pd.DataFrame({
        "spot_id": [f"spot_{i:04d}" for i in range(n_spots)],
        "hbam_score": rng.normal(0, 1, n_spots),
        "dysfunction_score": rng.uniform(0, 2, n_spots),
        "maturation_score": rng.uniform(0, 2, n_spots),
    })
    obs.index = obs["spot_id"]

    var = pd.DataFrame({"gene_names": gene_names}, index=gene_names)

    adata = ad.AnnData(
        X=X, obs=obs, var=var,
        obsm={"spatial": coords},
        layers={"raw": X.copy()},
    )
    return adata


@pytest.fixture
def sample_weights():
    """Create sample gene weights for testing."""
    return pd.DataFrame({
        "gene": ["ITIH2", "SELENBP1", "C3", "ACLY", "FASN",
                 "TXNRD1", "ETFA", "NDUFB4", "CD63", "MFAP4"],
        "weight": [0.1, 0.09, 0.08, 0.07, 0.06, 0.1, 0.09, 0.08, 0.07, 0.06],
        "category": ["dysfunction"] * 5 + ["maturation"] * 5,
        "modality": ["test"] * 10,
    })


def test_load_h5ad(tmp_path):
    """Test h5ad loader."""
    import anndata as ad
    from hbam.data.loaders import load_h5ad

    rng = np.random.default_rng(42)
    X = rng.normal(size=(10, 20)).astype(np.float32)
    coords = rng.uniform(0, 50, size=(10, 2)).astype(np.float32)

    adata = ad.AnnData(X=X, obsm={"spatial": coords})
    path = tmp_path / "test.h5ad"
    adata.write(path)

    loaded = load_h5ad(path)

    assert loaded.n_obs == 10
    assert loaded.n_vars == 20
    assert "spatial" in loaded.obsm
    assert "raw" in loaded.layers
    assert "gene_names" in loaded.var.columns


def test_compute_spatial_hbam_score(spatial_adata_with_scores, sample_weights):
    """Test per-spot HBAM scoring."""
    import anndata as ad
    from hbam.engine.score import compute_spatial_hbam_score

    # Remove pre-computed scores to test fresh computation
    adata = spatial_adata_with_scores.copy()
    del adata.obs["hbam_score"]
    del adata.obs["dysfunction_score"]
    del adata.obs["maturation_score"]

    result = compute_spatial_hbam_score(adata, sample_weights)

    assert "hbam_score" in result.obs.columns
    assert "dysfunction_score" in result.obs.columns
    assert "maturation_score" in result.obs.columns
    assert not result.obs["hbam_score"].isna().any()
    assert result.n_obs == adata.n_obs


def test_spatial_hbam_map(spatial_adata_with_scores, tmp_path):
    """Test spatial HBAM map figure generation."""
    from hbam.output.spatial_figures import fig_spatial_hbam_map

    paths = fig_spatial_hbam_map(spatial_adata_with_scores, tmp_path)

    assert len(paths) >= 1
    assert all(p.exists() for p in paths)


def test_spatial_gene_overlay(spatial_adata_with_scores, tmp_path):
    """Test spatial gene overlay figure."""
    from hbam.output.spatial_figures import fig_spatial_gene_overlay

    genes = ["ITIH2", "SELENBP1", "C3", "ACLY"]
    paths = fig_spatial_gene_overlay(spatial_adata_with_scores, genes, tmp_path)

    assert len(paths) >= 1
    assert all(p.exists() for p in paths)


def test_identify_region_biomarkers(spatial_adata_with_scores):
    """Test region-specific biomarker identification."""
    from hbam.translate.spatial_biomarkers import identify_region_biomarkers

    results = identify_region_biomarkers(spatial_adata_with_scores, top_n=10)

    assert isinstance(results, pd.DataFrame)
    assert "gene" in results.columns
    assert "log2fc" in results.columns
    assert "fdr" in results.columns
    assert "direction" in results.columns
    assert len(results) <= 10


def test_spatial_hbam_summary(spatial_adata_with_scores):
    """Test spatial HBAM summary statistics."""
    from hbam.translate.spatial_biomarkers import spatial_hbam_summary

    summary = spatial_hbam_summary(spatial_adata_with_scores)

    assert summary["n_spots"] == 200
    assert "mean" in summary
    assert "std" in summary
    assert summary["has_spatial"] is True


def test_spatial_no_coords(tmp_path):
    """Test graceful handling when spatial coords are missing."""
    import anndata as ad
    from hbam.output.spatial_figures import fig_spatial_hbam_map

    adata = ad.AnnData(X=np.random.default_rng(42).normal(size=(10, 5)).astype(np.float32))
    paths = fig_spatial_hbam_map(adata, tmp_path)

    assert len(paths) >= 1  # Should produce placeholder figure


def test_existing_tests_still_pass():
    """Regression check - verify core pipeline imports still work."""
    from hbam.pipeline import run_pipeline, PipelineResult
    from hbam.cli import app
    from hbam.config import load_config
    from hbam.data.loaders import load_maxquant, load_diann_matrix, load_h5ad
    from hbam.engine.score import compute_hbam_score, compute_spatial_hbam_score
    from hbam.output.spatial_figures import fig_spatial_hbam_map
    from hbam.translate.spatial_biomarkers import identify_region_biomarkers

    assert True  # All imports succeeded
