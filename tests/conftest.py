"""Shared test fixtures for the HBAM pipeline."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest


@pytest.fixture(scope="session")
def seed() -> int:
    """Fixed random seed for deterministic tests."""
    return 42


@pytest.fixture(scope="session")
def rng(seed: int) -> np.random.Generator:
    """Seeded random number generator."""
    return np.random.default_rng(seed)


@pytest.fixture(scope="session")
def make_proteomics_adata(rng: np.random.Generator):
    """Factory fixture for synthetic proteomics AnnData.

    Returns a function that creates AnnData with realistic proteomics structure:
    - Log-normal intensity distribution
    - ~10% missing values (MNAR-like, biased toward low intensity)
    - Sample metadata with condition and timepoint
    """
    def _make(
        n_samples: int = 20,
        n_proteins: int = 500,
        conditions: list[str] | None = None,
        seed: int = 42,
    ):
        import anndata as ad

        local_rng = np.random.default_rng(seed)

        if conditions is None:
            conditions = ["young", "old"]

        # Generate log-normal intensities (typical for proteomics)
        X = local_rng.lognormal(mean=20, sigma=2, size=(n_samples, n_proteins))

        # Inject condition effect in first 50 proteins (aging signature)
        n_de = min(50, n_proteins)
        condition_labels = [conditions[i % len(conditions)] for i in range(n_samples)]
        for i in range(n_samples):
            if condition_labels[i] == conditions[-1]:  # "old" condition
                X[i, :n_de] *= local_rng.lognormal(mean=0.5, sigma=0.3, size=n_de)

        # Inject MNAR missingness (~10%, biased toward low values)
        mask = local_rng.random(X.shape) < 0.1
        low_intensity = X < np.percentile(X, 30)
        mask = mask | (low_intensity & (local_rng.random(X.shape) < 0.15))
        X[mask] = np.nan

        # Create observation metadata
        obs = pd.DataFrame({
            "sample_id": [f"sample_{i:03d}" for i in range(n_samples)],
            "condition": condition_labels,
            "timepoint": [i // (n_samples // 4) + 1 for i in range(n_samples)],
            "batch": [f"batch_{i % 3}" for i in range(n_samples)],
        })
        obs.index = obs["sample_id"]

        # Create variable metadata
        gene_names = [f"GENE_{i:04d}" for i in range(n_proteins)]
        # Add some realistic gene names for the DE genes
        real_genes = [
            "TP53", "MTOR", "CDKN1A", "IL6", "TNF", "MYOD1", "PAX7", "MYH1",
            "TGFB1", "FOXO3", "SIRT1", "IGF1", "MSTN", "ATROGIN1", "MURF1",
        ]
        for i, g in enumerate(real_genes[:min(len(real_genes), n_de)]):
            gene_names[i] = g

        var = pd.DataFrame({
            "gene_names": gene_names,
            "protein_ids": [f"P{i:05d}" for i in range(n_proteins)],
        })
        var.index = gene_names

        adata = ad.AnnData(
            X=X.copy(),
            obs=obs,
            var=var,
            layers={"raw": X.copy()},
        )

        return adata

    return _make


@pytest.fixture(scope="session")
def make_spatial_adata(rng: np.random.Generator):
    """Factory fixture for synthetic spatial transcriptomics AnnData.

    Returns a function that creates AnnData mimicking Stereo-seq data:
    - Zero-inflated count distribution
    - Spatial coordinates in .obsm["spatial"]
    - Sparse expression patterns
    """
    def _make(
        n_bins: int = 200,
        n_genes: int = 300,
        seed: int = 42,
    ):
        import anndata as ad

        local_rng = np.random.default_rng(seed)

        # Generate zero-inflated negative binomial counts
        # (typical for spatial transcriptomics)
        mu = local_rng.gamma(shape=2, scale=3, size=(1, n_genes))
        X = local_rng.poisson(mu, size=(n_bins, n_genes)).astype(np.float32)

        # Zero inflation (~60% zeros, realistic for Stereo-seq)
        zero_mask = local_rng.random(X.shape) < 0.6
        X[zero_mask] = 0

        # Spatial coordinates (grid-like)
        grid_size = int(np.ceil(np.sqrt(n_bins)))
        coords = np.array(
            [(i % grid_size, i // grid_size) for i in range(n_bins)],
            dtype=np.float32,
        )

        # Spatially variable genes (first 20 genes have spatial pattern)
        n_sv = min(20, n_genes)
        for g in range(n_sv):
            center = local_rng.uniform(0, grid_size, size=2)
            dists = np.sqrt(((coords - center) ** 2).sum(axis=1))
            spatial_signal = np.exp(-dists / (grid_size / 3))
            X[:, g] = local_rng.poisson(mu[0, g] * (1 + 3 * spatial_signal))

        # Observation metadata
        obs = pd.DataFrame({
            "bin_id": [f"bin_{i:04d}" for i in range(n_bins)],
            "total_counts": X.sum(axis=1).astype(int),
            "n_genes_detected": (X > 0).sum(axis=1).astype(int),
            "region": [f"region_{i % 4}" for i in range(n_bins)],
        })
        obs.index = obs["bin_id"]

        # Variable metadata
        gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]
        spatial_genes = [
            "Ttn", "Myh1", "Myh2", "Acta1", "Ckm", "Mb", "Des", "Vim",
            "Col1a1", "Col3a1", "Fn1", "Lama2", "Pax7", "Myod1", "Myf5",
            "Sox9", "Acan", "Col2a1", "Runx2", "Sp7",
        ]
        for i, g in enumerate(spatial_genes[:min(len(spatial_genes), n_sv)]):
            gene_names[i] = g

        var = pd.DataFrame({
            "gene_names": gene_names,
        })
        var.index = gene_names

        adata = ad.AnnData(
            X=X,
            obs=obs,
            var=var,
            obsm={"spatial": coords},
            layers={"raw": X.copy()},
        )

        return adata

    return _make


@pytest.fixture(scope="session")
def make_mudata(make_proteomics_adata, make_spatial_adata):
    """Factory fixture for synthetic MuData with multiple modalities.

    Creates a MuData object with proteomics and spatial modalities
    that share some gene names (for integration testing).
    """
    def _make(
        n_samples: int = 20,
        n_proteins: int = 500,
        n_bins: int = 200,
        n_genes: int = 300,
        n_shared_genes: int = 150,
        seed: int = 42,
    ):
        import mudata as md

        prot = make_proteomics_adata(n_samples=n_samples, n_proteins=n_proteins, seed=seed)
        spat = make_spatial_adata(n_bins=n_bins, n_genes=n_genes, seed=seed + 1)

        # Ensure some genes are shared (for integration testing)
        shared_names = [f"SHARED_{i:04d}" for i in range(n_shared_genes)]
        prot_genes = list(prot.var_names)
        spat_genes = list(spat.var_names)

        # Replace some gene names with shared names
        n_to_replace = min(n_shared_genes, len(prot_genes), len(spat_genes))
        for i in range(n_to_replace):
            # Skip the first few (realistic names)
            prot_idx = min(15 + i, len(prot_genes) - 1)
            spat_idx = min(20 + i, len(spat_genes) - 1)
            prot_genes[prot_idx] = shared_names[i]
            spat_genes[spat_idx] = shared_names[i]

        prot.var_names = prot_genes
        prot.var["gene_names"] = prot_genes
        spat.var_names = spat_genes
        spat.var["gene_names"] = spat_genes

        mdata = md.MuData({"proteomics": prot, "spatial": spat})

        return mdata

    return _make


@pytest.fixture(scope="session")
def make_config():
    """Factory fixture for test PipelineConfig."""
    def _make(**overrides: Any):
        import sys
        sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
        from hbam.config import PipelineConfig

        # Sensible test defaults (smaller than production)
        defaults: dict[str, Any] = {
            "seed": 42,
            "log_level": "WARNING",  # Quiet during tests
            "integration": {"min_gene_overlap": 50, "n_factors": 5},
            "engine": {"bootstrap_iterations": 3, "bootstrap_fraction": 0.8},
            "translation": {"top_biomarkers": 10, "panel_sizes": [5, 10]},
            "output": {"figure_dpi": 72},  # Lower DPI for test speed
        }

        # Deep merge overrides
        def _merge(base: dict, updates: dict) -> dict:
            for k, v in updates.items():
                if isinstance(v, dict) and isinstance(base.get(k), dict):
                    base[k] = _merge(base[k], v)
                else:
                    base[k] = v
            return base

        merged = _merge(defaults, overrides)
        return PipelineConfig(**merged)

    return _make


@pytest.fixture
def sample_proteomics(make_proteomics_adata):
    """Pre-built proteomics AnnData for quick tests."""
    return make_proteomics_adata(n_samples=10, n_proteins=100, seed=42)


@pytest.fixture
def sample_spatial(make_spatial_adata):
    """Pre-built spatial AnnData for quick tests."""
    return make_spatial_adata(n_bins=50, n_genes=80, seed=42)


@pytest.fixture
def sample_mudata(make_mudata):
    """Pre-built MuData for quick tests."""
    return make_mudata(
        n_samples=10, n_proteins=100, n_bins=50, n_genes=80,
        n_shared_genes=30, seed=42,
    )


@pytest.fixture
def tmp_config(tmp_path, make_config):
    """Write a test config to a temp file and return path + config."""
    import yaml
    config = make_config()
    config_path = tmp_path / "test_config.yaml"
    data = config.model_dump(mode="json")
    with open(config_path, "w") as f:
        yaml.dump(data, f)
    return config_path, config


# ---------------------------------------------------------------------------
# Test data file fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def maxquant_file(tmp_path) -> Path:
    """Create a minimal MaxQuant proteinGroups.txt fixture."""
    data = {
        "Protein IDs": ["P12345", "P23456", "P34567", "REV__P99999", "CON__P88888"],
        "Gene names": ["TP53", "MTOR", "CDKN1A", "FAKE_REV", "FAKE_CON"],
        "Reverse": ["", "", "", "+", ""],
        "Potential contaminant": ["", "", "", "", "+"],
        "Only identified by site": ["", "", "", "", ""],
        "LFQ intensity S1": [1e6, 2e6, 3e6, 5e5, 1e5],
        "LFQ intensity S2": [1.1e6, 1.9e6, 3.2e6, 4e5, 2e5],
        "LFQ intensity S3": [0.9e6, 2.1e6, 2.8e6, 6e5, 1.5e5],
    }
    df = pd.DataFrame(data)
    filepath = tmp_path / "proteinGroups.txt"
    df.to_csv(filepath, sep="\t", index=False)
    return filepath


@pytest.fixture
def diann_file(tmp_path) -> Path:
    """Create a minimal DIA-NN report.tsv fixture."""
    rows = []
    rng = np.random.default_rng(42)
    for run in ["Run_01", "Run_02", "Run_03"]:
        for pg in ["P12345;P12346", "P23456", "P34567"]:
            rows.append({
                "Run": run,
                "Protein.Group": pg,
                "Genes": pg.split(";")[0].replace("P", "GENE"),
                "PG.MaxLFQ": rng.lognormal(20, 2),
            })
    df = pd.DataFrame(rows)
    filepath = tmp_path / "report.tsv"
    df.to_csv(filepath, sep="\t", index=False)
    return filepath


@pytest.fixture
def gem_file(tmp_path) -> Path:
    """Create a minimal Stereo-seq GEM fixture."""
    rng = np.random.default_rng(42)
    rows = []
    genes = ["Ttn", "Myh1", "Acta1", "Des", "Vim"]
    for gene in genes:
        n_spots = rng.integers(5, 20)
        for _ in range(n_spots):
            rows.append({
                "geneID": gene,
                "x": rng.integers(0, 100),
                "y": rng.integers(0, 100),
                "MIDCount": rng.integers(1, 50),
            })
    df = pd.DataFrame(rows)
    filepath = tmp_path / "stereo.gem"
    df.to_csv(filepath, sep="\t", index=False)
    return filepath
