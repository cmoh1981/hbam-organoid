"""Summary tables for the HBAM pipeline."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger

from hbam.utils.logging import log_step


def generate_summary_table(
    mudata,
    scores: pd.DataFrame,
    weights: pd.DataFrame,
) -> pd.DataFrame:
    """Generate pipeline summary statistics table.

    Args:
        mudata: Processed MuData.
        scores: HBAM scores.
        weights: Gene weights.

    Returns:
        Summary DataFrame.
    """
    with log_step("generate_summary_table"):
        rows = []

        # Modality info
        for name, mod in mudata.mod.items():
            rows.append({"Metric": f"Samples ({name})", "Value": str(mod.n_obs)})
            rows.append({"Metric": f"Genes ({name})", "Value": str(mod.n_vars)})

        # Score info
        rows.append({"Metric": "HBAM mean ± SD", "Value": f"{scores['hbam_score'].mean():.4f} ± {scores['hbam_score'].std():.4f}"})
        rows.append({"Metric": "HBAM range", "Value": f"[{scores['hbam_score'].min():.4f}, {scores['hbam_score'].max():.4f}]"})

        # Weight info
        rows.append({"Metric": "Total weighted genes", "Value": str(len(weights))})
        rows.append({"Metric": "Maturation genes", "Value": str((weights["category"] == "maturation").sum())})
        rows.append({"Metric": "Dysfunction genes", "Value": str((weights["category"] == "dysfunction").sum())})

        return pd.DataFrame(rows)


def generate_gene_table(weights: pd.DataFrame) -> pd.DataFrame:
    """Generate supplementary gene weight table.

    Args:
        weights: Gene weights DataFrame.

    Returns:
        Formatted table sorted by absolute weight.
    """
    df = weights.copy()
    df["abs_weight"] = df["weight"].abs()
    df = df.sort_values("abs_weight", ascending=False)
    df["rank"] = range(1, len(df) + 1)
    df = df.drop(columns=["abs_weight"])
    return df


def save_table(
    df: pd.DataFrame,
    name: str,
    output_dir: Path,
    formats: list[str] | None = None,
) -> list[Path]:
    """Save table in multiple formats.

    Args:
        df: DataFrame to save.
        name: Base filename.
        output_dir: Output directory.
        formats: List of formats (default: ["csv", "tex"]).

    Returns:
        List of saved paths.
    """
    if formats is None:
        formats = ["csv", "tex"]

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths = []
    for fmt in formats:
        path = output_dir / f"{name}.{fmt}"
        if fmt == "csv":
            df.to_csv(path, index=False)
        elif fmt == "tex":
            df.to_latex(path, index=False, escape=True)
        elif fmt == "tsv":
            df.to_csv(path, sep="\t", index=False)
        paths.append(path)
        logger.debug(f"Saved table: {path}")

    return paths
