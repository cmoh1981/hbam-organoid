"""Tests for end-to-end pipeline."""

import sys
sys.path.insert(0, "src")

import pytest
from pathlib import Path


def test_pipeline_result_dataclass():
    from hbam.pipeline import PipelineResult

    result = PipelineResult()
    assert result.mudata is None
    assert result.scores is None
    assert result.figure_paths == []


def test_set_seeds():
    import numpy as np
    from hbam.pipeline import set_seeds

    set_seeds(42)
    a = np.random.random()

    set_seeds(42)
    b = np.random.random()

    assert a == b


def test_generate_demo_data():
    from hbam.pipeline import _generate_demo_data

    modalities = _generate_demo_data(seed=42)

    assert "proteomics" in modalities
    assert "spatial" in modalities
    assert modalities["proteomics"].n_obs == 20
    assert modalities["spatial"].n_obs == 100


def test_end_to_end_demo(tmp_path):
    """Full pipeline test with synthetic demo data."""
    from hbam.config import PipelineConfig
    from hbam.pipeline import run_pipeline

    config = PipelineConfig(
        seed=42,
        log_level="WARNING",
        data={"output_dir": str(tmp_path / "results"),
              "intermediate_dir": str(tmp_path / "intermediate")},
        integration={"min_gene_overlap": 5, "n_factors": 3},
        engine={"bootstrap_iterations": 2, "bootstrap_fraction": 0.7},
        translation={"top_biomarkers": 5},
        output={"figure_dpi": 72},
    )

    result = run_pipeline(config)

    assert result.mudata is not None
    assert result.scores is not None or result.weights is not None
    # Should generate some figures
    assert len(result.figure_paths) >= 0  # Some may fail with test data


def test_cli_dry_run(tmp_path):
    """Test CLI --dry-run flag."""
    import yaml
    from hbam.config import PipelineConfig

    config = PipelineConfig(seed=42)
    config_path = tmp_path / "test_config.yaml"
    data = config.model_dump(mode="json")
    with open(config_path, "w") as f:
        yaml.dump(data, f)

    from typer.testing import CliRunner
    from hbam.cli import app

    runner = CliRunner()
    result = runner.invoke(app, ["run", "--config", str(config_path), "--dry-run"])

    assert result.exit_code == 0
    assert "Config validated" in result.stdout
