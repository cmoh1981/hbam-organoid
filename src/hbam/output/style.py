"""Publication-ready figure styling."""

from __future__ import annotations

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from loguru import logger

# Colorblind-safe palette (Okabe-Ito)
PALETTE = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "red": "#D55E00",
    "purple": "#CC79A7",
    "cyan": "#56B4E9",
    "yellow": "#F0E442",
    "black": "#000000",
    # Semantic mappings
    "young": "#0072B2",
    "old": "#D55E00",
    "aging": "#D55E00",
    "maturation": "#0072B2",
    "dysfunction": "#D55E00",
    "modality_1": "#0072B2",
    "modality_2": "#E69F00",
    "modality_3": "#009E73",
    "significant": "#D55E00",
    "non_significant": "#999999",
    "up": "#D55E00",
    "down": "#0072B2",
}

# Figure dimensions (inches)
SINGLE_COL = 3.5
DOUBLE_COL = 7.0
FULL_PAGE_HEIGHT = 9.0


def set_publication_style(
    font_size_label: int = 8,
    font_size_title: int = 10,
    font_family: str = "Arial",
) -> None:
    """Set matplotlib rcParams for publication-quality figures."""
    matplotlib.use("Agg")  # Non-interactive backend

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [font_family, "DejaVu Sans", "Helvetica"],
        "font.size": font_size_label,
        "axes.titlesize": font_size_title,
        "axes.labelsize": font_size_label,
        "xtick.labelsize": font_size_label - 1,
        "ytick.labelsize": font_size_label - 1,
        "legend.fontsize": font_size_label - 1,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "lines.linewidth": 1.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "pdf.fonttype": 42,  # TrueType fonts in PDF
        "ps.fonttype": 42,
    })


def save_figure(
    fig: plt.Figure,
    name: str,
    output_dir: Path,
    formats: list[str] | None = None,
    dpi: int = 300,
) -> list[Path]:
    """Save figure in multiple formats.

    Args:
        fig: Matplotlib figure.
        name: Base filename (without extension).
        output_dir: Output directory.
        formats: List of formats (default: ["pdf", "png"]).
        dpi: DPI for raster formats.

    Returns:
        List of saved file paths.
    """
    if formats is None:
        formats = ["pdf", "png"]

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths = []
    for fmt in formats:
        path = output_dir / f"{name}.{fmt}"
        fig.savefig(path, format=fmt, dpi=dpi, bbox_inches="tight")
        paths.append(path)
        logger.debug(f"Saved figure: {path}")

    plt.close(fig)
    return paths


def get_condition_colors(conditions: list[str]) -> dict[str, str]:
    """Map conditions to palette colors.

    Args:
        conditions: List of condition labels.

    Returns:
        Dict mapping condition -> color hex code.
    """
    color_list = [PALETTE["blue"], PALETTE["orange"], PALETTE["green"],
                  PALETTE["red"], PALETTE["purple"], PALETTE["cyan"]]

    colors = {}
    for i, cond in enumerate(sorted(conditions)):
        cond_lower = cond.lower()
        if cond_lower in PALETTE:
            colors[cond] = PALETTE[cond_lower]
        else:
            colors[cond] = color_list[i % len(color_list)]

    return colors
