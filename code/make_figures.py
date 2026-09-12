#!/usr/bin/env python3
"""Generate publication figures Fig 1-4 for the APAP zonation manuscript.

Follows the dataviz skill's validated categorical palette and mark specs:
thin 2px marks, recessive solid gridlines, ink-colored text, no dual axes,
legend for >=2 series with selective direct labels.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

sys.path.insert(0, str(Path(__file__).parent))
from zonated_apap_model import SCHEMES, POS, run, N_HEPATOCYTES  # noqa: E402
from intracellular_apap_model import dose_to_p0  # noqa: E402

FIGDIR = Path(__file__).resolve().parent.parent / "figures"

# --- validated categorical palette (light) + ink/chrome ---------------------
BLUE = "#2a78d6"
AQUA = "#1baf7a"
ORANGE = "#eb6834"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"

ASSUMED = ORANGE      # the replaced qualitative assumption
MEASURED = BLUE       # the measured recalibration
PERICENTRAL = BLUE
PERIPORTAL = AQUA

plt.rcParams.update({
    "font.size": 9,
    "axes.edgecolor": BASELINE,
    "axes.linewidth": 0.8,
    "axes.labelcolor": INK2,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.titleweight": "bold",
    "axes.titlesize": 9.5,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "legend.frameon": False,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

PARAM_LABEL = {
    "k450": "CYP450  (k450)",
    "kG": "Glucuronidation  (kG)",
    "kS": "Sulfation  (kS)",
    "kGSH": "GSH conjugation  (kGSH)",
    "bG": "GSH synthesis  (bG)",
}


def ramp(mp, mc):
    return mp + (mc - mp) * POS


def fig1_schematic():
    fig = plt.figure(figsize=(7.2, 2.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.7, 1], wspace=0.30)

    # --- left: sinusoid with 16 hepatocytes -------------------------------
    ax = fig.add_subplot(gs[0])
    ax.set_xlim(-0.6, 16.6)
    ax.set_ylim(0, 4.3)
    ax.axis("off")

    ax.add_patch(Rectangle((0, 3.4), 16, 0.42, facecolor=GRID, edgecolor=BASELINE, lw=0.8))
    ax.annotate("", xy=(16, 3.61), xytext=(0, 3.61),
                arrowprops=dict(arrowstyle="-|>", lw=1.6, color=MUTED))
    ax.text(8, 3.87, "sinusoidal blood flow", ha="center", va="bottom",
            fontsize=8, color=MUTED)

    for i in range(N_HEPATOCYTES):
        frac = i / (N_HEPATOCYTES - 1)
        ax.add_patch(Rectangle((i, 1.55), 0.92, 1.05,
                               facecolor=plt.cm.Blues(0.35 + 0.55 * frac),
                               edgecolor=BASELINE, lw=0.6))
    ax.text(0, 1.1, "1 (periportal)", ha="left", fontsize=7.5, color=MUTED)
    ax.text(16, 1.1, "16 (pericentral)", ha="right", fontsize=7.5, color=MUTED)

    ax.text(0, 2.85, "Zone 1\n(periportal)", ha="left", va="center", fontsize=7.5, color=INK2)
    ax.text(16, 2.85, "Zone 3\n(pericentral)", ha="right", va="center", fontsize=7.5, color=INK2)

    ax.text(0, 4.15, "portal vein", ha="left", fontsize=8, color=INK)
    ax.text(16, 4.15, "central vein", ha="right", fontsize=8, color=INK)

    # --- right: intracellular reaction tree --------------------------------
    ax = fig.add_subplot(gs[1])
    ax.set_xlim(0, 4.5)
    ax.set_ylim(0, 4.3)
    ax.axis("off")

    def node(x, y, s, fc="white", fs=7.5):
        ax.add_patch(FancyBboxPatch((x - 0.48, y - 0.27), 0.96, 0.54,
                                    boxstyle="round,pad=0.03",
                                    facecolor=fc, edgecolor=BASELINE, lw=0.7))
        ax.text(x, y, s, ha="center", va="center", fontsize=fs, color=INK)

    def arrow(x1, y1, x2, y2):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                     mutation_scale=7, lw=1.1, color=MUTED))

    def lab(x, y, s):
        ax.text(x, y, s, fontsize=6.5, color=MUTED, ha="center", va="center")

    node(2.25, 4.0, "APAP  (P)")
    node(0.8, 2.7, "APAP-sulfate")
    node(2.25, 2.7, "APAP-\nglucuronide")
    node(3.7, 2.7, "NAPQI  (N)")
    node(2.45, 1.15, "APAP-GSH\ndetoxified")
    node(3.9, 1.15, "adducts  (C)", fc="#f7d4d1")

    arrow(2.25, 3.73, 0.8, 2.97)
    arrow(2.25, 3.73, 2.25, 2.97)
    arrow(2.25, 3.73, 3.7, 2.97)
    arrow(3.7, 2.43, 2.45, 1.42)
    arrow(3.7, 2.43, 3.9, 1.42)

    lab(1.35, 3.5, "SULT")
    lab(2.25, 3.5, "UGT")
    lab(3.15, 3.5, "CYP450")
    lab(2.8, 1.95, "+ GSH")
    lab(4.15, 1.95, "kPSH")

    ax.set_title("Intracellular kinetics", fontsize=9, color=INK2, pad=6)

    fig.savefig(FIGDIR / "fig1_schematic.png")
    plt.close(fig)


def fig2_gradients():
    params = ["k450", "kG", "kS", "kGSH", "bG"]
    fig, axes = plt.subplots(1, 5, figsize=(10.5, 2.5), sharex=True, sharey=True)

    x = np.arange(1, N_HEPATOCYTES + 1)
    for ax, p in zip(axes, params):
        assumed = ramp(*SCHEMES["assumed"][p])
        measured = ramp(*SCHEMES["measured"][p])
        ax.plot(x, assumed, color=ASSUMED, ls="--", lw=2)
        ax.plot(x, measured, color=MEASURED, ls="-", lw=2)
        ax.axhline(1.0, color=BASELINE, lw=0.8)
        ax.set_yscale("log", base=2)
        ax.set_yticks([0.5, 1, 2, 4])
        ax.set_yticklabels(["0.5×", "1×", "2×", "4×"])
        ax.set_ylim(0.25, 5)
        ax.set_title(PARAM_LABEL[p], fontsize=8)
        ax.grid(True, which="both", lw=0.6)
        ax.set_xticks([1, 16])
        ax.set_xticklabels(["portal", "central"], fontsize=7.5)
        ax.tick_params(axis="x", pad=1)
        ax.tick_params(axis="y", labelsize=7.5)

    axes[0].set_ylabel("relative activity\n(fold vs uniform)", fontsize=8)

    # highlight the flip on kGSH
    axes[3].annotate("direction\nreversed", xy=(16, 1.85), xytext=(7, 3.0),
                     fontsize=7, color=INK2, ha="center",
                     arrowprops=dict(arrowstyle="->", lw=0.8, color=MUTED))

    # shared legend via proxy
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=ASSUMED, ls="--", lw=2, label="assumed (Gebhardt-derived)"),
               Line2D([], [], color=MEASURED, ls="-", lw=2, label="measured (2026 human)")]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, -0.02))

    fig.savefig(FIGDIR / "fig2_gradients.png")
    plt.close(fig)


def fig3_damage():
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9), sharey=True)
    x = np.arange(1, N_HEPATOCYTES + 1)

    for ax, dose in zip(axes, (4.0, 16.0)):
        p0 = dose_to_p0(dose)
        for scheme, color, ls, lab in [("assumed", ASSUMED, "--", "assumed"),
                                        ("measured", MEASURED, "-", "measured")]:
            t, Y, _ = run(scheme, dose)
            Cpk = Y[:, :, 4].max(axis=0)
            ax.plot(x, Cpk / p0 * 100, color=color, ls=ls, lw=2)
        ax.set_title(f"{dose:.0f} g", fontsize=9.5)
        ax.set_xlabel("hepatocyte (1 = periportal, 16 = pericentral)", fontsize=8)
        ax.set_xticks([1, 4, 8, 12, 16])
        ax.tick_params(labelsize=8)
        ax.grid(True, lw=0.6)

    axes[0].set_ylabel("peak adduct fraction\nof dose (%)", fontsize=8)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=ASSUMED, ls="--", lw=2, label="assumed"),
               Line2D([], [], color=MEASURED, ls="-", lw=2, label="measured")]
    fig.legend(handles=handles, loc="upper left", fontsize=8, frameon=False,
               bbox_to_anchor=(0.12, 0.98))

    fig.savefig(FIGDIR / "fig3_damage.png")
    plt.close(fig)


def fig4_uncertainty():
    npz_path = FIGDIR.parent / "results" / "uncertainty_changes.npz"
    if not npz_path.exists():
        print("  (skip fig4: run code/uncertainty.py first)")
        return
    data = np.load(npz_path)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8))

    for ax, dose, key in zip(axes, (4.0, 16.0), ("changes_4g", "changes_16g")):
        changes = data[key]
        med = np.median(changes)
        lo, hi = np.percentile(changes, [5, 95])
        ax.hist(changes, bins=26, color=BLUE, alpha=0.75, edgecolor="white", lw=0.3)
        ax.axvline(0, color=INK2, lw=1)
        ax.axvline(med, color=ORANGE, lw=1.5, ls="--")
        ax.set_title(f"{dose:.0f} g", fontsize=9.5)
        ax.set_xlabel("Δ pericentral adducts (%)", fontsize=8)
        ax.tick_params(labelsize=8)
        ax.grid(True, lw=0.6, axis="y")
        prot = (changes < 0).mean() * 100
        ax.text(0.97, 0.96, f"median {med:+.0f}%\n90% interval [{lo:+.0f}, {hi:+.0f}]\n"
                            f"{prot:.0f}% protective",
                transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=INK2)

    axes[0].set_ylabel("Monte Carlo draws", fontsize=8)
    fig.savefig(FIGDIR / "fig4_uncertainty.png")
    plt.close(fig)


def main():
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig1_schematic()
    fig2_gradients()
    fig3_damage()
    fig4_uncertainty()
    for f in sorted(FIGDIR.glob("*.png")):
        print("wrote", f.resolve())


if __name__ == "__main__":
    main()
