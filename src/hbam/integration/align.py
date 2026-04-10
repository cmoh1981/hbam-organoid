"""Gene alignment and feature scaling across modalities."""

from __future__ import annotations

import mudata as md
import numpy as np
import pandas as pd
from loguru import logger

from hbam.utils.logging import log_filter, log_metric, log_step
from hbam.utils.validation import ValidationGate, validate_gene_overlap


def align_genes(
    mudata: md.MuData,
    min_overlap: int = 1000,
) -> md.MuData:
    """Subset all modalities to shared genes.

    Args:
        mudata: Input MuData with multiple modalities.
        min_overlap: Minimum required gene overlap.

    Returns:
        MuData with all modalities having identical var_names.

    Raises:
        PipelineValidationError: If overlap < min_overlap.
    """
    with log_step("align_genes", n_modalities=len(mudata.mod)):
        # Compute gene sets per modality
        gene_sets = {}
        for name, mod in mudata.mod.items():
            genes = set(mod.var_names.tolist())
            gene_sets[name] = genes
            log_metric(f"genes_{name}", len(genes))

        # Compute intersection
        all_genes = list(gene_sets.values())
        shared = all_genes[0]
        for gs in all_genes[1:]:
            shared = shared & gs

        shared_sorted = sorted(shared)
        n_shared = len(shared_sorted)

        log_metric("shared_genes", n_shared)

        # Validate overlap
        for name_a, genes_a in gene_sets.items():
            for name_b, genes_b in gene_sets.items():
                if name_a < name_b:
                    gate = ValidationGate(f"gene_overlap_{name_a}_{name_b}")
                    gate.add(validate_gene_overlap(genes_a, genes_b, min_overlap=min_overlap))
                    gate.enforce()

        # Subset each modality
        result_mods = {}
        for name, mod in mudata.mod.items():
            before = mod.n_vars
            subset = mod[:, mod.var_names.isin(shared_sorted)].copy()
            # Reorder to match
            subset = subset[:, shared_sorted].copy()
            result_mods[name] = subset
            log_filter(f"align_{name}", before=before, after=subset.n_vars, reason="gene alignment")

        result = md.MuData(result_mods)

        return result


def scale_features(
    mudata: md.MuData,
    method: str = "zscore",
) -> md.MuData:
    """Standardize features across modalities.

    Args:
        mudata: Input MuData (aligned genes).
        method: "zscore" (mean=0, std=1) or "minmax" (0-1 range).

    Returns:
        Scaled MuData with .layers["unscaled"] preserved.
    """
    with log_step("scale_features", method=method):
        for name, mod in mudata.mod.items():
            X = mod.X.copy().astype(np.float64)

            # Store unscaled
            mod.layers["unscaled"] = X.copy()

            if method == "zscore":
                means = np.nanmean(X, axis=0)
                stds = np.nanstd(X, axis=0)
                stds[stds < 1e-10] = 1.0  # Avoid division by zero
                X = (X - means) / stds
            elif method == "minmax":
                mins = np.nanmin(X, axis=0)
                maxs = np.nanmax(X, axis=0)
                ranges = maxs - mins
                ranges[ranges < 1e-10] = 1.0
                X = (X - mins) / ranges

            mod.X = X
            log_metric(f"scaled_{name}_mean", f"{np.nanmean(X):.4f}")
            log_metric(f"scaled_{name}_std", f"{np.nanstd(X):.4f}")

        return mudata


def report_alignment(mudata: md.MuData) -> dict:
    """Generate alignment statistics report.

    Args:
        mudata: Aligned MuData.

    Returns:
        Dict with alignment statistics.
    """
    gene_sets = {name: set(mod.var_names) for name, mod in mudata.mod.items()}

    all_genes = set()
    for gs in gene_sets.values():
        all_genes |= gs

    shared = set.intersection(*gene_sets.values()) if gene_sets else set()

    report = {
        "n_modalities": len(mudata.mod),
        "total_unique_genes": len(all_genes),
        "shared_genes": len(shared),
        "per_modality": {
            name: {"n_genes": len(gs), "n_samples": mudata.mod[name].n_obs}
            for name, gs in gene_sets.items()
        },
    }

    # Jaccard similarities
    names = list(gene_sets.keys())
    for i, n1 in enumerate(names):
        for j, n2 in enumerate(names):
            if i < j:
                union = len(gene_sets[n1] | gene_sets[n2])
                inter = len(gene_sets[n1] & gene_sets[n2])
                jaccard = inter / union if union > 0 else 0
                report[f"jaccard_{n1}_{n2}"] = jaccard

    return report
