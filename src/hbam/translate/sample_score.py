"""Sample-level HBAM scoring and condition comparison."""

from __future__ import annotations

import mudata as md
import numpy as np
import pandas as pd
from loguru import logger
from scipy.stats import mannwhitneyu

from hbam.utils.logging import log_metric, log_step


def score_samples(
    mudata: md.MuData,
    weights: pd.DataFrame,
) -> pd.DataFrame:
    """Apply HBAM weights to compute per-sample scores.

    Args:
        mudata: Input MuData.
        weights: Gene weights DataFrame.

    Returns:
        DataFrame with sample_id, hbam_score, dysfunction_score, maturation_score,
        plus any metadata columns from .obs.
    """
    with log_step("score_samples"):
        first_mod = list(mudata.mod.keys())[0]
        mod = mudata.mod[first_mod]

        gene_to_idx = {g: i for i, g in enumerate(mod.var_names)}

        dysfunction = np.zeros(mod.n_obs)
        maturation = np.zeros(mod.n_obs)

        for _, row in weights.iterrows():
            gene = row["gene"]
            if gene not in gene_to_idx:
                gene = next((g for g in gene_to_idx if g.upper() == gene.upper()), None)
                if gene is None:
                    continue

            idx = gene_to_idx[gene]
            vals = np.nan_to_num(mod.X[:, idx].copy(), nan=0.0)

            if row["category"] == "dysfunction":
                dysfunction += vals * row["weight"]
            else:
                maturation += vals * row["weight"]

        result = pd.DataFrame({
            "sample_id": mod.obs_names,
            "hbam_score": dysfunction - maturation,
            "dysfunction_score": dysfunction,
            "maturation_score": maturation,
        })

        # Add metadata
        for col in mod.obs.columns:
            if col not in result.columns:
                result[col] = mod.obs[col].values

        result = result.set_index("sample_id")

        log_metric("scored_samples", len(result))
        log_metric("score_range", f"[{result['hbam_score'].min():.4f}, {result['hbam_score'].max():.4f}]")

        return result


def compare_conditions(
    scores: pd.DataFrame,
    group_col: str = "condition",
) -> dict:
    """Statistical comparison of HBAM scores between groups.

    Args:
        scores: DataFrame with hbam_score and group column.
        group_col: Column defining groups.

    Returns:
        Dict with p_value, effect_size, group_means, group_medians.
    """
    with log_step("compare_conditions", group_col=group_col):
        if group_col not in scores.columns:
            logger.warning(f"Column '{group_col}' not found")
            return {"p_value": None, "effect_size": None}

        groups = scores[group_col].unique()
        if len(groups) < 2:
            return {"p_value": None, "effect_size": None, "n_groups": len(groups)}

        group_data = {g: scores[scores[group_col] == g]["hbam_score"].values for g in groups}

        # Pairwise comparison (first two groups)
        g1, g2 = sorted(groups)[:2]
        d1, d2 = group_data[g1], group_data[g2]

        stat, p = mannwhitneyu(d1, d2, alternative="two-sided")

        # Cohen's d
        pooled_std = np.sqrt((np.var(d1) + np.var(d2)) / 2)
        cohens_d = (np.mean(d2) - np.mean(d1)) / pooled_std if pooled_std > 0 else 0

        result = {
            "p_value": float(p),
            "statistic": float(stat),
            "effect_size_cohens_d": float(cohens_d),
            "group_means": {str(g): float(np.mean(v)) for g, v in group_data.items()},
            "group_medians": {str(g): float(np.median(v)) for g, v in group_data.items()},
            "group_sizes": {str(g): len(v) for g, v in group_data.items()},
            "comparison": f"{g2}_vs_{g1}",
        }

        log_metric("comparison_pvalue", f"{p:.6f}")
        log_metric("comparison_cohens_d", f"{cohens_d:.4f}")

        return result
