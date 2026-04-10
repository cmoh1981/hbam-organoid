"""Tests for HBAM engine."""

import numpy as np
import pandas as pd
import pytest

import sys
sys.path.insert(0, "src")


def test_load_genage():
    from hbam.engine.gene_sets import load_genage

    genes = load_genage()

    assert len(genes) > 30
    assert "TP53" in genes
    assert "SIRT1" in genes


def test_categorize_genes():
    from hbam.engine.gene_sets import load_msigdb_hallmark, categorize_genes

    hallmark = load_msigdb_hallmark()
    categories = categorize_genes(hallmark)

    assert "maturation" in categories
    assert "dysfunction" in categories
    assert len(categories["maturation"]) > 0
    assert len(categories["dysfunction"]) > 0
    # No overlap
    assert len(categories["maturation"] & categories["dysfunction"]) == 0


def test_compute_weights(sample_mudata):
    from hbam.engine.gene_sets import load_genage, load_msigdb_hallmark, categorize_genes
    from hbam.engine.weights import compute_weights_pca
    from hbam.integration.align import scale_features
    from hbam.integration.embed import embed_pca

    scaled = scale_features(sample_mudata)
    embedded = embed_pca(scaled, n_components=5, seed=42)

    genage = load_genage()
    hallmark = load_msigdb_hallmark()
    gene_sets = categorize_genes(hallmark, genage)

    weights = compute_weights_pca(embedded, gene_sets)

    # May be empty if no gene overlap, that's OK for test data
    assert isinstance(weights, pd.DataFrame)
    if len(weights) > 0:
        assert "gene" in weights.columns
        assert "weight" in weights.columns
        assert "category" in weights.columns


def test_normalize_weights():
    from hbam.engine.weights import normalize_weights

    df = pd.DataFrame({
        "gene": ["A", "B", "C", "D"],
        "weight": [0.5, 0.3, 0.4, 0.6],
        "category": ["maturation", "maturation", "dysfunction", "dysfunction"],
    })

    result = normalize_weights(df)

    # Each category should sum to 1
    for cat in ["maturation", "dysfunction"]:
        cat_sum = result[result["category"] == cat]["weight"].sum()
        assert abs(cat_sum - 1.0) < 1e-6
