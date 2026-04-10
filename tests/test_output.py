"""Tests for output layer."""

import numpy as np
import pytest
from pathlib import Path

import sys
sys.path.insert(0, "src")


def test_style_system():
    from hbam.output.style import set_publication_style, PALETTE, get_condition_colors

    set_publication_style()

    assert "blue" in PALETTE
    assert "dysfunction" in PALETTE
    assert "maturation" in PALETTE

    colors = get_condition_colors(["young", "old"])
    assert "young" in colors
    assert "old" in colors


def test_fig_pipeline_overview(tmp_path):
    from hbam.output.figures import fig_pipeline_overview

    paths = fig_pipeline_overview(tmp_path)

    assert len(paths) >= 1
    assert all(p.exists() for p in paths)


def test_fig_hbam_distribution(tmp_path):
    import pandas as pd
    from hbam.output.figures import fig_hbam_distribution

    scores = pd.DataFrame({
        "hbam_score": np.random.default_rng(42).normal(0, 1, 20),
        "condition": ["young"] * 10 + ["old"] * 10,
    })

    paths = fig_hbam_distribution(scores, tmp_path)
    assert len(paths) >= 1
    assert all(p.exists() for p in paths)


def test_fig_weight_heatmap(tmp_path):
    import pandas as pd
    from hbam.output.figures import fig_weight_heatmap

    weights = pd.DataFrame({
        "gene": [f"GENE_{i}" for i in range(20)],
        "weight": np.random.default_rng(42).uniform(-1, 1, 20),
        "category": ["maturation"] * 10 + ["dysfunction"] * 10,
    })

    paths = fig_weight_heatmap(weights, tmp_path, top_n=10)
    assert len(paths) >= 1
    assert all(p.exists() for p in paths)


def test_tables():
    import pandas as pd
    from hbam.output.tables import generate_gene_table

    weights = pd.DataFrame({
        "gene": ["A", "B", "C"],
        "weight": [0.5, 0.3, 0.2],
        "category": ["maturation", "dysfunction", "maturation"],
    })

    result = generate_gene_table(weights)
    assert "rank" in result.columns
    assert result.iloc[0]["weight"] == 0.5  # Highest weight first


def test_save_table(tmp_path):
    import pandas as pd
    from hbam.output.tables import save_table

    df = pd.DataFrame({"col1": [1, 2], "col2": ["a", "b"]})
    paths = save_table(df, "test", tmp_path, formats=["csv"])

    assert len(paths) == 1
    assert paths[0].exists()
    assert paths[0].suffix == ".csv"
