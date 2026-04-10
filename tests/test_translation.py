"""Tests for translation layer."""

import numpy as np
import pandas as pd
import pytest
import sys
sys.path.insert(0, "src")


def test_score_samples(sample_mudata):
    from hbam.translate.sample_score import score_samples
    from hbam.integration.align import scale_features
    from hbam.integration.embed import embed_pca
    from hbam.engine.gene_sets import load_genage, load_msigdb_hallmark, categorize_genes
    from hbam.engine.weights import compute_weights_pca, normalize_weights

    scaled = scale_features(sample_mudata)
    embedded = embed_pca(scaled, n_components=5, seed=42)

    genage = load_genage()
    hallmark = load_msigdb_hallmark()
    gene_sets = categorize_genes(hallmark, genage)

    weights = compute_weights_pca(embedded, gene_sets)

    if len(weights) > 0:
        weights = normalize_weights(weights)
        scores = score_samples(embedded, weights)

        assert "hbam_score" in scores.columns
        assert "dysfunction_score" in scores.columns
        assert "maturation_score" in scores.columns
        assert len(scores) > 0
        assert not scores["hbam_score"].isna().any()


def test_compare_conditions():
    from hbam.translate.sample_score import compare_conditions

    scores = pd.DataFrame(
        {
            "hbam_score": np.concatenate([
                np.random.default_rng(42).normal(0, 1, 10),
                np.random.default_rng(42).normal(2, 1, 10),
            ]),
            "condition": ["young"] * 10 + ["old"] * 10,
        }
    )

    result = compare_conditions(scores, group_col="condition")

    assert "p_value" in result
    assert "effect_size_cohens_d" in result
    assert result["p_value"] is not None
    assert result["p_value"] < 0.05  # Should be significant with this effect size


def test_select_biomarkers():
    from hbam.translate.biomarkers import select_biomarkers

    weights = pd.DataFrame({
        "gene": [f"GENE_{i}" for i in range(50)],
        "weight": np.random.default_rng(42).uniform(0, 1, 50),
        "category": ["maturation"] * 25 + ["dysfunction"] * 25,
    })

    panel = select_biomarkers(weights, n=10)

    assert len(panel) == 10
    assert "rank" in panel.columns
    # Sorted by absolute weight descending
    abs_weights = panel["weight"].abs().values
    assert abs_weights[0] >= abs_weights[-1]


def test_compute_effect_sizes():
    from hbam.translate.clinical import compute_effect_sizes

    scores = pd.DataFrame({
        "hbam_score": np.concatenate([
            np.random.default_rng(42).normal(0, 0.5, 15),
            np.random.default_rng(43).normal(3, 0.5, 15),
        ]),
        "condition": ["young"] * 15 + ["old"] * 15,
    })

    result = compute_effect_sizes(scores)

    assert len(result) > 0
    assert "cohens_d" in result.columns
    assert "hedges_g" in result.columns
    assert "p_value" in result.columns
    assert abs(result.iloc[0]["cohens_d"]) > 1.0  # Large effect


def test_discrimination_analysis():
    from hbam.translate.clinical import discrimination_analysis

    scores = pd.DataFrame({
        "hbam_score": np.concatenate([
            np.random.default_rng(42).normal(0, 1, 20),
            np.random.default_rng(42).normal(3, 1, 20),
        ]),
        "condition": ["young"] * 20 + ["old"] * 20,
    })

    result = discrimination_analysis(scores, positive_label="old")

    assert "auc" in result
    assert result["auc"] > 0.7  # Should discriminate well
    assert "sensitivity" in result
    assert "specificity" in result


def test_clinical_summary():
    from hbam.translate.clinical import generate_clinical_summary

    scores = pd.DataFrame({
        "hbam_score": [1.0, 2.0, 3.0],
        "condition": ["young", "old", "old"],
    })
    effect = pd.DataFrame({
        "comparison": ["old_vs_young"],
        "cohens_d": [1.5],
        "hedges_g": [1.4],
        "p_value": [0.01],
        "n1": [1],
        "n2": [2],
        "mean_1": [1.0],
        "mean_2": [2.5],
    })
    disc = {
        "auc": 0.85,
        "sensitivity": 0.8,
        "specificity": 0.9,
        "optimal_threshold": 1.5,
    }

    summary = generate_clinical_summary(scores, effect, disc)

    assert isinstance(summary, str)
    assert "HBAM" in summary
    assert "ROC" in summary or "AUC" in summary
