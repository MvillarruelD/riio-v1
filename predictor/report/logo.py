"""Shared sequence-logo rendering and report colours for the live HTML report."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D

# Editable text in SVG exports, with broadly available fallbacks.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.8,
                     "axes.spines.right": False, "axes.spines.top": False,
                     "legend.frameon": False})

NMI = {"baseline_dark": "#484878", "baseline_mid": "#7884B4", "baseline_soft": "#B4C0E4",
       "ours_base": "#E4CCD8", "ours_large": "#F0C0CC", "bg_lilac": "#E0E0F0",
       "bg_aqua": "#E0F0F0", "bg_peach": "#F0E0D0", "neutral_light": "#D8D8D8",
       "neutral_mid": "#A8A8A8", "neutral_dark": "#606060", "delta_up": "#2E9E44",
       "delta_down": "#E53935"}
BASE_COLORS = {"A": "#2E9E44", "C": "#3775BA", "G": "#E2A12C", "T": "#D9544D"}
_MONO = FontProperties(family="monospace", weight="bold")


def draw_logo(ax, pwm, importance, *, ymax=2.0):
    """Draw an information-content logo on ``ax`` from a 4 x L PWM."""
    pwm = np.asarray(pwm)
    importance = np.asarray(importance)
    length = pwm.shape[1]
    for position in range(length):
        order = np.argsort(pwm[:, position])
        y = 0.0
        for index in order:
            base = "ACGT"[index]
            height = float(pwm[index, position]) * float(importance[position])
            if height <= 0.01:
                y += height
                continue
            text_path = TextPath((0, 0), base, size=1, prop=_MONO)
            bounds = text_path.get_extents()
            transform = (
                Affine2D()
                .translate(-bounds.x0, -bounds.y0)
                .scale(1.0 / max(bounds.width, 1e-6), height / max(bounds.height, 1e-6))
                .translate(position + 0.06, y)
                + ax.transData
            )
            ax.add_patch(PathPatch(text_path, transform=transform, color=BASE_COLORS[base], lw=0))
            y += height
    ax.set_xlim(0, length)
    ax.set_ylim(0, ymax)
    ax.set_ylabel("bits", fontsize=7.5)
    # Ticks on column centres, numbered from 1 like the sequence itself.
    ticks = sorted({1, (length + 1) // 2, length})
    ax.set_xticks([t - 1 + 0.56 for t in ticks])
    ax.set_xticklabels([str(t) for t in ticks])
    ax.set_xlabel("position in motif (bp)", fontsize=7.5)
    ax.tick_params(labelsize=7)
