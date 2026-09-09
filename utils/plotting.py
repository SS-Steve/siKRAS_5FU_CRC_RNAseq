"""Unified plotting style for the siKRAS RNA-seq pipeline."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ── Unified colour scheme for all figures ──────────────────────────────────
COLOR_UP = "#B31A21"      # upregulated / high group (red)
COLOR_DOWN = "#1465AC"    # downregulated / low group (blue)
COLOR_NS = "#A9A9A9"      # not significant (grey)
COLOR_HK = "#FFFFFF"      # housekeeping fill
COLOR_HK_EDGE = "#000000"
COLOR_GROUP_PALETTE = "Set3"

# Blue-white-red diverging colormap for heatmaps
BWR_CMAP = LinearSegmentedColormap.from_list(
    "BWR", [COLOR_DOWN, "#FFFFFF", COLOR_UP]
)


def setup_style() -> None:
    """Apply a consistent matplotlib style."""
    plt.rcParams.update({
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.titleweight": "bold",
        "axes.labelsize": 13,
        "axes.labelweight": "bold",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 150,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
    })


def save_figure(fig, path: str) -> None:
    """Save a figure to path, creating parent directories."""
    from pathlib import Path
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def style_legend(ax, **kwargs) -> None:
    """Apply a consistent legend style (white background + border)."""
    defaults = dict(fontsize=10, handletextpad=0.5, frameon=True,
                    facecolor='white', edgecolor='black')
    defaults.update(kwargs)
    ax.legend(**defaults)


def add_grid(ax, axis="both", alpha=0.3) -> None:
    """Add a light grey grid behind the data."""
    ax.grid(color="gray", alpha=alpha, axis=axis)
    ax.set_axisbelow(True)
