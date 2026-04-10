"""HBAM composite score computation and validation."""

from __future__ import annotations

from dataclasses import dataclass, field

import mudata as md
import numpy as np
import pandas as pd
from loguru import logger
from scipy.stats import spearmanr, mannwhitneyu

from hbam.utils.logging import log_metric, log_step


@dataclass
class BootstrapResult:
    """Result of bootstrap validation."""
    correlations: list[float]
    mean_correlation: float
    passed: bool
    threshold: float
    n_iterations: int


def compute_hbam_score(
    mudata: md.MuData,
    weights: pd.DataFrame,
) -> pd.Series:
    """Compute HBAM composite score per sample.

    HBAM = dysfunction_score - maturation_score
    Higher HBAM = more aging burden.

    Args:
        mudata: Input MuData.
        weights: DataFrame with gene, weight, category columns.

    Returns:
        Series of HBAM scores indexed by sample.
    """
    with log_step("compute_hbam_score"):
        # Use first modality for scoring (or concatenated)
        first_mod_name = list(mudata.mod.keys())[0]
        mod = mudata.mod[first_mod_name]

        gene_to_idx = {g: i for i, g in enumerate(mod.var_names)}

        dysfunction_score = np.zeros(mod.n_obs)
        maturation_score = np.zeros(mod.n_obs)

        n_dys_found = 0
        n_mat_found = 0

        for _, row in weights.iterrows():
            gene = row["gene"]
            weight = row["weight"]
            category = row["category"]

            if gene not in gene_to_idx:
                # Try uppercase
                gene_upper = gene.upper()
                match = next((g for g in gene_to_idx if g.upper() == gene_upper), None)
                if match:
                    gene = match
                else:
                    continue

            idx = gene_to_idx[gene]
            values = mod.X[:, idx].copy()

            # Replace NaN with 0 for scoring
            if np.issubdtype(values.dtype, np.floating):
                values = np.nan_to_num(values, nan=0.0)

            if category == "dysfunction":
                dysfunction_score += values * weight
                n_dys_found += 1
            elif category == "maturation":
                maturation_score += values * weight
                n_mat_found += 1

        hbam = dysfunction_score - maturation_score

        # Store in MuData
        mod.obs["hbam_score"] = hbam
        mod.obs["dysfunction_score"] = dysfunction_score
        mod.obs["maturation_score"] = maturation_score

        log_metric("hbam_mean", f"{np.mean(hbam):.4f}")
        log_metric("hbam_std", f"{np.std(hbam):.4f}")
        log_metric("genes_used_dysfunction", n_dys_found)
        log_metric("genes_used_maturation", n_mat_found)

        scores = pd.Series(hbam, index=mod.obs_names, name="hbam_score")
        return scores


def validate_score_direction(
    scores: pd.Series,
    metadata: pd.DataFrame,
    age_col: str = "condition",
    old_label: str = "old",
) -> bool:
    """Validate that HBAM positively correlates with aging.

    Args:
        scores: HBAM scores.
        metadata: Sample metadata.
        age_col: Column with age/condition info.
        old_label: Label for the aged condition.

    Returns:
        True if direction is correct (higher HBAM = more aging).
    """
    if age_col not in metadata.columns:
        logger.warning(f"Column '{age_col}' not in metadata. Cannot validate direction.")
        return True

    conditions = metadata[age_col].unique()
    if len(conditions) < 2:
        return True

    # Check if "old" samples have higher HBAM
    old_mask = metadata[age_col] == old_label
    if not old_mask.any():
        return True

    old_scores = scores[old_mask].mean()
    young_scores = scores[~old_mask].mean()

    correct = old_scores > young_scores

    if not correct:
        logger.warning(
            f"HBAM direction may be inverted: old={old_scores:.4f}, young={young_scores:.4f}. "
            "Consider flipping sign."
        )
    else:
        logger.info(f"HBAM direction validated: old={old_scores:.4f} > young={young_scores:.4f}")

    return correct


def bootstrap_validate(
    mudata: md.MuData,
    gene_sets: dict[str, set[str]],
    compute_weights_fn,
    n_iter: int = 5,
    fraction: float = 0.8,
    threshold: float = 0.8,
    seed: int = 42,
) -> BootstrapResult:
    """Validate HBAM weight stability via bootstrap resampling.

    Args:
        mudata: Input MuData.
        gene_sets: Gene set categories.
        compute_weights_fn: Function to compute weights (mofa or pca).
        n_iter: Number of bootstrap iterations.
        fraction: Fraction of samples per iteration.
        threshold: Minimum Spearman correlation required.
        seed: Random seed.

    Returns:
        BootstrapResult with correlations and pass/fail.
    """
    with log_step("bootstrap_validate", n_iter=n_iter, fraction=fraction):
        # Compute full-data weights
        full_weights = compute_weights_fn(mudata, gene_sets)

        if len(full_weights) == 0:
            logger.warning("No weights computed. Bootstrap validation skipped.")
            return BootstrapResult(
                correlations=[], mean_correlation=0.0,
                passed=False, threshold=threshold, n_iterations=0,
            )

        full_weight_vec = full_weights.set_index("gene")["weight"]

        correlations = []
        first_mod = list(mudata.mod.keys())[0]
        n_samples = mudata.mod[first_mod].n_obs
        n_subsample = max(3, int(n_samples * fraction))

        for i in range(n_iter):
            rng = np.random.default_rng(seed + i)
            indices = rng.choice(n_samples, size=n_subsample, replace=False)

            # Subset mudata
            sub_mods = {}
            for name, mod in mudata.mod.items():
                if indices.max() < mod.n_obs:
                    sub_mods[name] = mod[indices].copy()
                else:
                    sub_idx = indices[indices < mod.n_obs]
                    sub_mods[name] = mod[sub_idx].copy()

            import mudata as md_lib
            sub_mudata = md_lib.MuData(sub_mods)

            # Recompute weights
            try:
                boot_weights = compute_weights_fn(sub_mudata, gene_sets)
                if len(boot_weights) == 0:
                    continue

                boot_weight_vec = boot_weights.set_index("gene")["weight"]

                # Compute correlation on shared genes
                shared = full_weight_vec.index.intersection(boot_weight_vec.index)
                if len(shared) < 3:
                    continue

                r, _ = spearmanr(full_weight_vec[shared], boot_weight_vec[shared])
                if np.isfinite(r):
                    correlations.append(float(r))
            except Exception as e:
                logger.warning(f"Bootstrap iteration {i} failed: {e}")
                continue

        mean_corr = np.mean(correlations) if correlations else 0.0
        passed = all(c >= threshold for c in correlations) if correlations else False

        result = BootstrapResult(
            correlations=correlations,
            mean_correlation=float(mean_corr),
            passed=passed,
            threshold=threshold,
            n_iterations=len(correlations),
        )

        log_metric("bootstrap_mean_spearman", f"{mean_corr:.4f}")
        log_metric("bootstrap_passed", passed)

        return result


def bootstrap_confidence_intervals(
    mudata: md.MuData,
    gene_sets: dict[str, set[str]],
    compute_weights_fn,
    n_iter: int = 100,
    fraction: float = 0.8,
    seed: int = 42,
    ci: float = 0.95,
) -> pd.DataFrame:
    """Compute bootstrap confidence intervals for gene weights.

    Args:
        mudata: Input MuData.
        gene_sets: Gene set categories.
        compute_weights_fn: Function to compute weights.
        n_iter: Number of bootstrap iterations.
        fraction: Fraction of samples per iteration.
        seed: Random seed.
        ci: Confidence interval level (default 0.95).

    Returns:
        DataFrame with gene, weight_mean, weight_lower, weight_upper, weight_std.
    """
    with log_step("bootstrap_confidence_intervals", n_iter=n_iter):
        all_weights: dict[str, list[float]] = {}

        first_mod = list(mudata.mod.keys())[0]
        n_samples = mudata.mod[first_mod].n_obs
        n_subsample = max(3, int(n_samples * fraction))

        for i in range(n_iter):
            rng = np.random.default_rng(seed + i)
            indices = rng.choice(n_samples, size=n_subsample, replace=False)

            sub_mods = {}
            for name, mod in mudata.mod.items():
                sub_idx = indices[indices < mod.n_obs]
                if len(sub_idx) >= 3:
                    sub_mods[name] = mod[sub_idx].copy()

            if not sub_mods:
                continue

            import mudata as md_lib
            sub_mudata = md_lib.MuData(sub_mods)

            try:
                boot_weights = compute_weights_fn(sub_mudata, gene_sets)
                for _, row in boot_weights.iterrows():
                    gene = row["gene"]
                    if gene not in all_weights:
                        all_weights[gene] = []
                    all_weights[gene].append(row["weight"])
            except Exception:
                continue

        if not all_weights:
            return pd.DataFrame(columns=["gene", "weight_mean", "weight_lower", "weight_upper", "weight_std"])

        alpha = (1 - ci) / 2
        rows = []
        for gene, weights_list in all_weights.items():
            arr = np.array(weights_list)
            rows.append({
                "gene": gene,
                "weight_mean": float(np.mean(arr)),
                "weight_lower": float(np.percentile(arr, alpha * 100)),
                "weight_upper": float(np.percentile(arr, (1 - alpha) * 100)),
                "weight_std": float(np.std(arr)),
                "n_bootstraps": len(arr),
            })

        result = pd.DataFrame(rows).sort_values("weight_mean", ascending=False)
        log_metric("bootstrap_ci_genes", len(result))

        return result
