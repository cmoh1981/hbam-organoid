"""Tests for QC, normalization, and imputation."""

import numpy as np
import pytest

import sys
sys.path.insert(0, "src")


def test_filter_missingness(sample_proteomics):
    from hbam.data.qc import filter_missingness

    result = filter_missingness(sample_proteomics, threshold=0.5)

    # Should have fewer or equal features
    assert result.n_vars <= sample_proteomics.n_vars
    assert result.n_obs == sample_proteomics.n_obs


def test_filter_variance(sample_proteomics):
    from hbam.data.qc import filter_variance

    result = filter_variance(sample_proteomics, percentile=10)

    assert result.n_vars < sample_proteomics.n_vars
    assert result.n_obs == sample_proteomics.n_obs


def test_qc_report(sample_proteomics):
    from hbam.data.qc import qc_report

    report = qc_report(sample_proteomics)

    assert "n_samples" in report
    assert "n_features" in report
    assert "missing_fraction" in report
    assert report["n_samples"] == sample_proteomics.n_obs
    assert report["n_features"] == sample_proteomics.n_vars


def test_normalize_proteomics(sample_proteomics):
    from hbam.data.normalize import normalize_proteomics

    result = normalize_proteomics(sample_proteomics)

    assert "raw" in result.layers
    # Should be log-transformed (much smaller values)
    assert np.nanmean(result.X) < np.nanmean(result.layers["raw"])


def test_normalize_spatial(sample_spatial):
    from hbam.data.normalize import normalize_spatial

    result = normalize_spatial(sample_spatial, method="total")

    assert "raw" in result.layers
    assert result.X is not None


def test_impute_auto(sample_proteomics):
    from hbam.data.impute import impute_auto

    result = impute_auto(sample_proteomics, seed=42)

    assert "imputed_mask" in result.layers
    # No NaN after imputation
    assert not np.isnan(result.X).any()


def test_detect_missingness_type(sample_proteomics):
    from hbam.data.impute import detect_missingness_type

    result = detect_missingness_type(sample_proteomics)

    assert result in ("MNAR", "MCAR", "MIXED")
