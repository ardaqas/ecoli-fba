"""Shared matplotlib style so every figure in the project looks the same."""

from pathlib import Path

import matplotlib as mpl
from matplotlib.figure import Figure

from ecoli_fba import RESULTS_DIR

FIGURES_DIR = RESULTS_DIR / "figures"

# Categorical colours, always used in this order (colour-blind-safe as
# neighbours). Text stays in the ink colours, never in a series colour.
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
NEUTRAL = "#c3c2b7"  # de-emphasised marks

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"


def apply_style() -> None:
    """Set matplotlib defaults: light surface, thin marks, recessive grid."""
    mpl.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "figure.dpi": 110,
            "savefig.dpi": 200,
            "savefig.bbox": "tight",
            "font.family": "sans-serif",
            "font.size": 10,
            "text.color": INK,
            "axes.titlesize": 11,
            "axes.titleweight": "bold",
            "axes.titlelocation": "left",
            "axes.labelcolor": INK_SECONDARY,
            "axes.edgecolor": AXIS,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "axes.axisbelow": True,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "grid.linestyle": "-",
            "xtick.color": INK_MUTED,
            "ytick.color": INK_MUTED,
            "xtick.labelcolor": INK_SECONDARY,
            "ytick.labelcolor": INK_SECONDARY,
            "legend.frameon": False,
            "patch.edgecolor": SURFACE,  # thin surface gap between bars
            "patch.linewidth": 1.0,
            "patch.force_edgecolor": True,
            "lines.linewidth": 2.0,
        }
    )


def save_figure(fig: Figure, name: str) -> Path:
    """Save a figure as PNG under results/figures/ and return its path."""
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    path = FIGURES_DIR / f"{name}.png"
    fig.savefig(path)
    return path
