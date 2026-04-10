"""Functional analysis for muscle proteomics data."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
from loguru import logger
from scipy import stats
from statsmodels.stats.multitest import multipletests

from hbam.utils.logging import log_filter, log_metric, log_step


def run_functional_analysis(
    adata: ad.AnnData,
    condition_col: str = "condition",
    reference: str | None = None,
    resolve_groups: bool = True,
) -> ad.AnnData:
    """Run functional annotation and differential expression on muscle proteomics.

    Args:
        adata: Input AnnData.
        condition_col: Column in .obs for DE comparison.
        reference: Reference condition for fold-change calculation.
        resolve_groups: Whether to resolve multi-gene protein groups.

    Returns:
        AnnData with results in .var:
        - de_pvalue, de_fdr, log2fc: differential expression results
        - functional_category: broad functional category (if annotatable)
    """
    with log_step("functional_analysis", condition_col=condition_col):
        result = adata.copy()

        # Resolve protein groups (semicolon-separated gene names)
        if resolve_groups and "protein_ids" in result.var.columns:
            result = _resolve_protein_groups(result)

        # Differential expression
        if condition_col in result.obs.columns:
            conditions = result.obs[condition_col].unique()

            if len(conditions) >= 2:
                if reference is None:
                    reference = sorted(conditions)[0]

                test_conditions = [c for c in conditions if c != reference]

                # Use first non-reference condition for primary comparison
                test_cond = test_conditions[0]

                ref_mask = result.obs[condition_col] == reference
                test_mask = result.obs[condition_col] == test_cond

                pvalues = np.ones(result.n_vars)
                log2fc = np.zeros(result.n_vars)

                X = result.X

                for j in range(result.n_vars):
                    ref_vals = X[ref_mask, j]
                    test_vals = X[test_mask, j]

                    # Remove NaN
                    ref_clean = ref_vals[~np.isnan(ref_vals)] if np.issubdtype(ref_vals.dtype, np.floating) else ref_vals
                    test_clean = test_vals[~np.isnan(test_vals)] if np.issubdtype(test_vals.dtype, np.floating) else test_vals

                    if len(ref_clean) >= 2 and len(test_clean) >= 2:
                        try:
                            stat, p = stats.mannwhitneyu(ref_clean, test_clean, alternative="two-sided")
                            pvalues[j] = p
                        except ValueError:
                            pass

                        ref_mean = np.mean(ref_clean)
                        test_mean = np.mean(test_clean)

                        if ref_mean > 0:
                            log2fc[j] = np.log2(test_mean / ref_mean) if test_mean > 0 else 0.0

                # FDR correction
                reject, fdr, _, _ = multipletests(pvalues, method="fdr_bh")

                result.var["de_pvalue"] = pvalues
                result.var["de_fdr"] = fdr
                result.var["log2fc"] = log2fc
                result.var["de_comparison"] = f"{test_cond}_vs_{reference}"

                n_sig = (fdr < 0.05).sum()
                n_up = ((fdr < 0.05) & (log2fc > 0)).sum()
                n_down = ((fdr < 0.05) & (log2fc < 0)).sum()

                log_metric("de_significant", int(n_sig), comparison=f"{test_cond}_vs_{reference}")
                log_metric("de_up", int(n_up))
                log_metric("de_down", int(n_down))
            else:
                logger.warning(f"Only one condition found in '{condition_col}', skipping DE")
                result.var["de_pvalue"] = 1.0
                result.var["de_fdr"] = 1.0
                result.var["log2fc"] = 0.0
        else:
            logger.warning(f"Condition column '{condition_col}' not found, skipping DE")
            result.var["de_pvalue"] = 1.0
            result.var["de_fdr"] = 1.0
            result.var["log2fc"] = 0.0

        # Simple functional categorization based on gene name patterns
        result.var["functional_category"] = _categorize_genes(result.var_names.tolist())

        return result


def _resolve_protein_groups(adata: ad.AnnData) -> ad.AnnData:
    """Resolve multi-gene protein groups by selecting first gene.

    When protein IDs contain semicolons (e.g., "P12345;P23456"),
    the gene name may also have multiple entries. We take the first.
    """
    gene_names = adata.var_names.tolist()
    new_names = []
    for name in gene_names:
        if ";" in name:
            new_names.append(name.split(";")[0].strip())
        else:
            new_names.append(name)

    n_resolved = sum(1 for old, new in zip(gene_names, new_names) if old != new)
    if n_resolved > 0:
        log_filter("resolve_protein_groups", before=len(gene_names), after=len(gene_names), reason=f"resolved {n_resolved} multi-gene groups")

    # Make unique
    seen = {}
    unique_names = []
    for name in new_names:
        if name in seen:
            seen[name] += 1
            unique_names.append(f"{name}_{seen[name]}")
        else:
            seen[name] = 0
            unique_names.append(name)

    adata = adata.copy()
    adata.var_names = unique_names
    adata.var["gene_names"] = unique_names

    return adata


def _categorize_genes(gene_names: list[str]) -> list[str]:
    """Simple functional categorization based on gene name patterns."""
    categories = {
        "muscle": ["MYH", "MYL", "ACTA", "TTN", "DES", "MYOD", "MYF", "PAX7", "CKM", "MB"],
        "inflammatory": ["IL", "TNF", "NFKB", "CXCL", "CCL", "STAT3", "JAK"],
        "metabolic": ["AMPK", "MTOR", "AKT", "IGF", "INS", "GLUT", "PPARG"],
        "senescence": ["CDKN", "TP53", "RB1", "TERT", "SIRT"],
        "ecm": ["COL", "FN1", "LAMA", "MMP", "TIMP"],
    }

    result = []
    for gene in gene_names:
        gene_upper = gene.upper()
        assigned = "other"
        for cat, prefixes in categories.items():
            if any(gene_upper.startswith(p) for p in prefixes):
                assigned = cat
                break
        result.append(assigned)

    return result
