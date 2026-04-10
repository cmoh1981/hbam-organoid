"""Clinical metrics and effect size computation."""

from __future__ import annotations

import numpy as np
import pandas as pd
from loguru import logger
from scipy.stats import mannwhitneyu

from hbam.utils.logging import log_metric, log_step


def compute_effect_sizes(
    scores: pd.DataFrame,
    comparisons: list[tuple[str, str]] | None = None,
    group_col: str = "condition",
    score_col: str = "hbam_score",
) -> pd.DataFrame:
    """Compute effect sizes for pairwise condition comparisons.

    Args:
        scores: DataFrame with scores and group column.
        comparisons: List of (group1, group2) pairs. None = all pairs.
        group_col: Column with group labels.
        score_col: Column with scores.

    Returns:
        DataFrame with comparison, cohens_d, hedges_g, p_value, n1, n2.
    """
    with log_step("compute_effect_sizes"):
        groups = scores[group_col].unique()

        if comparisons is None:
            comparisons = [(g1, g2) for i, g1 in enumerate(sorted(groups))
                          for g2 in sorted(groups)[i+1:]]

        rows = []
        for g1, g2 in comparisons:
            d1 = scores[scores[group_col] == g1][score_col].values
            d2 = scores[scores[group_col] == g2][score_col].values

            if len(d1) < 2 or len(d2) < 2:
                continue

            # Cohen's d
            pooled_std = np.sqrt(((len(d1)-1)*np.var(d1, ddof=1) + (len(d2)-1)*np.var(d2, ddof=1))
                                / (len(d1) + len(d2) - 2))
            cohens_d = (np.mean(d2) - np.mean(d1)) / pooled_std if pooled_std > 0 else 0

            # Hedge's g (bias-corrected)
            correction = 1 - 3 / (4 * (len(d1) + len(d2)) - 9)
            hedges_g = cohens_d * correction

            # Mann-Whitney U test
            stat, p = mannwhitneyu(d1, d2, alternative="two-sided")

            rows.append({
                "comparison": f"{g2}_vs_{g1}",
                "cohens_d": float(cohens_d),
                "hedges_g": float(hedges_g),
                "p_value": float(p),
                "n1": len(d1),
                "n2": len(d2),
                "mean_1": float(np.mean(d1)),
                "mean_2": float(np.mean(d2)),
            })

        return pd.DataFrame(rows)


def discrimination_analysis(
    scores: pd.DataFrame,
    group_col: str = "condition",
    positive_label: str | None = None,
    score_col: str = "hbam_score",
) -> dict:
    """ROC-AUC analysis for discrimination between conditions.

    Args:
        scores: DataFrame with scores and group column.
        group_col: Column with group labels.
        positive_label: Which group is the positive class.
        score_col: Column with scores.

    Returns:
        Dict with auc, optimal_threshold, sensitivity, specificity.
    """
    with log_step("discrimination_analysis"):
        from sklearn.metrics import roc_auc_score, roc_curve

        groups = sorted(scores[group_col].unique())
        if len(groups) < 2:
            return {"auc": None, "error": "Need at least 2 groups"}

        if positive_label is None:
            positive_label = groups[-1]

        binary = (scores[group_col] == positive_label).astype(int).values
        score_vals = scores[score_col].values

        auc = roc_auc_score(binary, score_vals)
        fpr, tpr, thresholds = roc_curve(binary, score_vals)

        # Youden's J
        j_idx = np.argmax(tpr - fpr)

        result = {
            "auc": float(auc),
            "optimal_threshold": float(thresholds[j_idx]),
            "sensitivity": float(tpr[j_idx]),
            "specificity": float(1 - fpr[j_idx]),
            "positive_label": positive_label,
            "n_positive": int(binary.sum()),
            "n_negative": int((1 - binary).sum()),
        }

        log_metric("auc", f"{auc:.4f}")
        log_metric("sensitivity", f"{result['sensitivity']:.4f}")
        log_metric("specificity", f"{result['specificity']:.4f}")

        return result


def generate_clinical_summary(
    scores: pd.DataFrame,
    effect_sizes: pd.DataFrame,
    discrimination: dict,
) -> str:
    """Generate formatted clinical summary text.

    Args:
        scores: Sample scores DataFrame.
        effect_sizes: Effect size results.
        discrimination: Discrimination analysis results.

    Returns:
        Formatted summary string for manuscript.
    """
    lines = ["## Clinical Summary", ""]

    # Score distribution
    hbam = scores["hbam_score"]
    lines.append(f"HBAM scores ranged from {hbam.min():.3f} to {hbam.max():.3f} "
                 f"(mean: {hbam.mean():.3f} ± {hbam.std():.3f}).")
    lines.append("")

    # Effect sizes
    if len(effect_sizes) > 0:
        for _, row in effect_sizes.iterrows():
            lines.append(
                f"**{row['comparison']}**: Cohen's d = {row['cohens_d']:.3f}, "
                f"Hedge's g = {row['hedges_g']:.3f}, "
                f"p = {row['p_value']:.2e} "
                f"(n₁={row['n1']}, n₂={row['n2']})"
            )
        lines.append("")

    # Discrimination
    if discrimination.get("auc"):
        lines.append(
            f"ROC-AUC: {discrimination['auc']:.3f} | "
            f"Sensitivity: {discrimination['sensitivity']:.1%} | "
            f"Specificity: {discrimination['specificity']:.1%} "
            f"(optimal threshold: {discrimination['optimal_threshold']:.3f})"
        )

    return "\n".join(lines)
