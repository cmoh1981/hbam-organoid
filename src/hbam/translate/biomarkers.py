"""Biomarker panel selection and evaluation."""

from __future__ import annotations

import mudata as md
import numpy as np
import pandas as pd
from loguru import logger

from hbam.utils.logging import log_metric, log_step


def select_biomarkers(
    weights: pd.DataFrame,
    n: int = 20,
) -> pd.DataFrame:
    """Select top N genes as biomarker panel by absolute weight.

    Args:
        weights: Gene weights DataFrame.
        n: Number of biomarkers to select.

    Returns:
        DataFrame with gene, weight, rank, category.
    """
    with log_step("select_biomarkers", n=n):
        df = weights.copy()
        df["abs_weight"] = df["weight"].abs()
        df = df.sort_values("abs_weight", ascending=False).head(n)
        df["rank"] = range(1, len(df) + 1)
        df = df.drop(columns=["abs_weight"])

        n_mat = (df["category"] == "maturation").sum()
        n_dys = (df["category"] == "dysfunction").sum()

        log_metric("biomarker_panel_size", len(df))
        log_metric("biomarker_maturation", int(n_mat))
        log_metric("biomarker_dysfunction", int(n_dys))

        return df


def evaluate_panel(
    mudata: md.MuData,
    panel_genes: list[str],
    weights: pd.DataFrame,
    group_col: str = "condition",
) -> dict:
    """Evaluate discrimination power of a biomarker panel.

    Args:
        mudata: Input MuData.
        panel_genes: List of panel gene names.
        weights: Full gene weights for comparison.
        group_col: Condition column.

    Returns:
        Dict with auc, sensitivity, specificity, panel_correlation.
    """
    with log_step("evaluate_panel", n_genes=len(panel_genes)):
        from hbam.translate.sample_score import score_samples

        # Score using only panel genes
        panel_weights = weights[weights["gene"].isin(panel_genes)].copy()
        panel_scores = score_samples(mudata, panel_weights)

        # Score using all genes for correlation
        full_scores = score_samples(mudata, weights)

        # Correlation between panel and full scores
        from scipy.stats import spearmanr
        shared_idx = panel_scores.index.intersection(full_scores.index)
        r, _ = spearmanr(panel_scores.loc[shared_idx, "hbam_score"],
                         full_scores.loc[shared_idx, "hbam_score"])

        result = {"panel_correlation": float(r), "n_genes": len(panel_genes)}

        # AUC if we have conditions
        first_mod = list(mudata.mod.keys())[0]
        mod = mudata.mod[first_mod]

        if group_col in mod.obs.columns:
            conditions = mod.obs[group_col].unique()
            if len(conditions) >= 2:
                try:
                    from sklearn.metrics import roc_auc_score

                    sorted_conds = sorted(conditions)
                    binary = (mod.obs[group_col] == sorted_conds[-1]).astype(int).values
                    scores_arr = panel_scores["hbam_score"].values[:len(binary)]

                    auc = roc_auc_score(binary, scores_arr)
                    result["auc"] = float(auc)

                    # Optimal threshold (Youden's J)
                    from sklearn.metrics import roc_curve
                    fpr, tpr, thresholds = roc_curve(binary, scores_arr)
                    j_idx = np.argmax(tpr - fpr)

                    result["optimal_threshold"] = float(thresholds[j_idx])
                    result["sensitivity"] = float(tpr[j_idx])
                    result["specificity"] = float(1 - fpr[j_idx])

                except Exception as e:
                    logger.warning(f"AUC calculation failed: {e}")

        return result


def optimize_panel_size(
    mudata: md.MuData,
    weights: pd.DataFrame,
    sizes: list[int] | None = None,
    group_col: str = "condition",
) -> pd.DataFrame:
    """Evaluate panel performance across different sizes.

    Args:
        mudata: Input MuData.
        weights: Full gene weights.
        sizes: Panel sizes to evaluate.
        group_col: Condition column.

    Returns:
        DataFrame with size, correlation, auc for each panel size.
    """
    if sizes is None:
        sizes = [5, 10, 15, 20, 30, 50]

    with log_step("optimize_panel_size", sizes=sizes):
        rows = []

        for n in sizes:
            panel = select_biomarkers(weights, n=n)
            panel_genes = panel["gene"].tolist()

            metrics = evaluate_panel(mudata, panel_genes, weights, group_col)

            rows.append({
                "panel_size": n,
                "correlation": metrics.get("panel_correlation", np.nan),
                "auc": metrics.get("auc", np.nan),
                "sensitivity": metrics.get("sensitivity", np.nan),
                "specificity": metrics.get("specificity", np.nan),
            })

        result = pd.DataFrame(rows)
        return result
