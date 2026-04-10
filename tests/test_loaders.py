"""Tests for data loaders."""

import numpy as np
import pytest

import sys
sys.path.insert(0, "src")


def test_load_maxquant(maxquant_file):
    from hbam.data.loaders import load_maxquant

    adata = load_maxquant(maxquant_file)

    # Should filter contaminants and reverse hits
    assert adata.n_obs == 3  # 3 samples (S1, S2, S3)
    assert adata.n_vars == 3  # 5 proteins - 1 reverse - 1 contaminant = 3
    assert "gene_names" in adata.var.columns
    assert "raw" in adata.layers


def test_load_diann(diann_file):
    from hbam.data.loaders import load_diann

    adata = load_diann(diann_file)

    assert adata.n_obs == 3  # 3 runs
    assert adata.n_vars == 3  # 3 protein groups
    assert "raw" in adata.layers


def test_load_matrix(tmp_path):
    import pandas as pd
    from hbam.data.loaders import load_matrix

    # Create test CSV
    data = {"gene": ["A", "B", "C"], "S1": [1.0, 2.0, 3.0], "S2": [4.0, 5.0, 6.0]}
    df = pd.DataFrame(data)
    path = tmp_path / "test.csv"
    df.to_csv(path, index=False)

    adata = load_matrix(path)
    assert adata.n_vars > 0
    assert "raw" in adata.layers


def test_load_stereo_gem(gem_file):
    from hbam.data.loaders import load_stereo_gem

    adata = load_stereo_gem(gem_file)

    assert adata.n_obs > 0
    assert adata.n_vars > 0
    assert "spatial" in adata.obsm
    assert adata.obsm["spatial"].shape[1] == 2  # x, y coordinates
    assert "raw" in adata.layers
