"""Edge case tests for pipeline robustness."""

import numpy as np
import pandas as pd
import pytest
import sys
sys.path.insert(0, "src")


def test_high_missingness_rejection():
    """Pipeline should filter features with high missingness."""
    import anndata as ad
    from hbam.data.qc import filter_missingness

    rng = np.random.default_rng(42)
    X = rng.normal(size=(10, 50)).astype(np.float64)
    X[rng.random(X.shape) < 0.6] = np.nan  # ~60% missing

    obs = pd.DataFrame({"sample_id": [f"S{i}" for i in range(10)]})
    obs.index = obs["sample_id"]
    var = pd.DataFrame({"gene_names": [f"G{i}" for i in range(50)]})
    var.index = var["gene_names"]

    adata = ad.AnnData(X=X, obs=obs, var=var)

    # threshold=0.9 means keep only features present in >=90% of samples.
    # With ~60% missing, most features have ~40% presence, so most are filtered.
    result = filter_missingness(adata, threshold=0.9)
    assert result.n_vars < adata.n_vars


def test_single_modality_pipeline():
    """Pipeline components should work with a single modality."""
    import mudata as md
    import anndata as ad
    from hbam.integration.align import scale_features
    from hbam.integration.embed import embed_pca

    rng = np.random.default_rng(42)
    X = rng.normal(size=(15, 30))
    obs = pd.DataFrame({"condition": ["young"] * 7 + ["old"] * 8})
    obs.index = [f"S{i}" for i in range(15)]
    var = pd.DataFrame(index=[f"G{i}" for i in range(30)])

    adata = ad.AnnData(X=X, obs=obs, var=var)
    mudata = md.MuData({"single": adata})

    scaled = scale_features(mudata)
    result = embed_pca(scaled, n_components=3, seed=42)

    assert "pca_factors" in result.obsm
    assert result.obsm["pca_factors"].shape == (15, 3)


def test_no_gene_overlap_error():
    """align_genes should raise PipelineValidationError when overlap is below min."""
    import mudata as md
    import anndata as ad
    from hbam.integration.align import align_genes
    from hbam.utils.validation import PipelineValidationError

    rng = np.random.default_rng(42)

    a1 = ad.AnnData(
        X=rng.normal(size=(5, 10)),
        var=pd.DataFrame(index=[f"GENE_A_{i}" for i in range(10)]),
    )
    a2 = ad.AnnData(
        X=rng.normal(size=(5, 10)),
        var=pd.DataFrame(index=[f"GENE_B_{i}" for i in range(10)]),
    )

    mudata = md.MuData({"mod1": a1, "mod2": a2})

    # min_overlap=5 but there is 0 overlap — should raise
    with pytest.raises(PipelineValidationError):
        align_genes(mudata, min_overlap=5)


def test_validation_gate_failure():
    """ValidationGate.enforce() should raise PipelineValidationError on failed checks."""
    from hbam.utils.validation import ValidationGate, validate_gene_overlap, PipelineValidationError

    gate = ValidationGate("test")
    gate.add(validate_gene_overlap({"A", "B"}, {"C", "D"}, min_overlap=3))

    assert not gate.passed

    with pytest.raises(PipelineValidationError):
        gate.enforce()


def test_imputation_preserves_mask():
    """impute_mnar should track which values were imputed in layers['imputed_mask']."""
    import anndata as ad
    from hbam.data.impute import impute_mnar

    rng = np.random.default_rng(42)
    X = rng.normal(10, 2, size=(10, 20))
    X[0, 0] = np.nan
    X[3, 5] = np.nan
    X[7, 15] = np.nan

    obs = pd.DataFrame(index=[f"S{i}" for i in range(10)])
    var = pd.DataFrame(index=[f"G{i}" for i in range(20)])
    adata = ad.AnnData(X=X, obs=obs, var=var)

    result = impute_mnar(adata, seed=42)

    assert "imputed_mask" in result.layers
    mask = result.layers["imputed_mask"]
    assert mask[0, 0] == 1.0   # Was NaN, now imputed
    assert mask[3, 5] == 1.0
    assert mask[7, 15] == 1.0
    assert mask[1, 1] == 0.0   # Was not NaN
    assert not np.isnan(result.X).any()


def test_config_invalid_values():
    """PipelineConfig should reject invalid parameter values via Pydantic."""
    from hbam.config import PipelineConfig, QCConfig
    from pydantic import ValidationError

    # missingness_threshold has ge=0.0 constraint — negative value must fail
    with pytest.raises(ValidationError):
        QCConfig(missingness_threshold=-0.5)

    # log_level validator rejects unknown strings
    with pytest.raises(ValidationError):
        PipelineConfig(log_level="INVALID")


def test_deterministic_pipeline():
    """_generate_demo_data with the same seed must produce identical results."""
    from hbam.pipeline import _generate_demo_data

    data1 = _generate_demo_data(seed=42)
    data2 = _generate_demo_data(seed=42)

    np.testing.assert_array_equal(data1["proteomics"].X, data2["proteomics"].X)
    np.testing.assert_array_equal(data1["spatial"].X, data2["spatial"].X)
