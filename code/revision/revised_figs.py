"""Revised Fig. 2 (readable, vector+raster) and new Fig. S1 (mechanism check for kG)."""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.integrate import solve_ivp

REPO = Path(sys.argv[1]) / "code"; OUT = Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REPO))
import intracellular_apap_model as im
import zonated_apap_model as zm

plt.rcParams.update({"font.size": 11, "axes.labelsize": 12, "axes.titlesize": 12,
                     "legend.fontsize": 11, "xtick.labelsize": 10.5, "ytick.labelsize": 10.5,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})
OR, BL, GR = "#eb6834", "#2a78d6", "#1baf7a"
x = np.arange(1, 17)

# ---------------- Fig 2 ----------------
labels = {"k450": "CYP450 ($k_{450}$)", "kG": "Glucuronidation ($k_G$)", "kS": "Sulfation ($k_S$)",
          "kGSH": "GSH conjugation ($k_{GSH}$)", "bG": "GSH synthesis ($b_G$)"}
folds = {"k450": ("5", "3.4"), "kG": ("9", "1.4"), "kS": ("1/9", "1.0"), "kGSH": ("1/9", "1.2"), "bG": ("1/9", "1.0")}
fig, axes = plt.subplots(2, 3, figsize=(10, 6.2), sharex=True, sharey=True)
for ax, p in zip(axes.flat, labels):
    a = zm.ramp(*zm.SCHEMES["assumed"][p]); m = zm.ramp(*zm.SCHEMES["measured"][p])
    ax.plot(x, a, color=OR, ls="--", lw=2.4); ax.plot(x, m, color=BL, lw=2.4)
    ax.axhline(1, color="0.6", lw=1)
    ax.set_yscale("log", base=2); ax.set_ylim(0.1, 2.2)
    ax.set_yticks([0.125, 0.25, 0.5, 1, 2]); ax.set_yticklabels(["0.125", "0.25", "0.5", "1", "2"])
    ax.set_title(labels[p]); ax.grid(alpha=0.3)
    ax.text(0.04, 0.05, f"assumed {folds[p][0]}$\\times$\nmeasured {folds[p][1]}$\\times$",
            transform=ax.transAxes, fontsize=9.5, va="bottom")
    ax.set_xticks([1, 8, 16])
for ax in axes[1]: ax.set_xlabel("hepatocyte (1 = portal, 16 = central)", fontsize=10)
for ax in axes[:, 0]: ax.set_ylabel("multiplier on baseline rate")
leg = axes.flat[5]; leg.axis("off")
leg.legend(handles=[Line2D([], [], color=OR, ls="--", lw=2.4, label="Assumed (Gebhardt-derived, prior models)"),
                    Line2D([], [], color=BL, lw=2.4, label="Measured (Weiss et al. 2026 bin profiles)")],
           loc="center", frameon=False)
fig.tight_layout()
for ext in ("pdf", "svg"): fig.savefig(OUT / f"fig2_gradients_revised.{ext}")
fig.savefig(OUT / "fig2_gradients_revised.png", dpi=600); plt.close(fig)

# ---------------- Fig S1: mechanism of the kG effect ----------------
def sim(scheme, dose):
    p = zm.build_parameters(scheme); p0 = im.dose_to_p0(dose)
    y0 = np.stack([np.full(16, p0), np.full(16, im.B_S / im.D_S), np.zeros(16), p["bG"] / im.D_G, np.zeros(16)], 1).ravel()
    t = np.concatenate([[0], np.geomspace(1e-5, 3.0, 4000)])
    s = solve_ivp(zm.make_rhs(p), (0, 3), y0, method="Radau", t_eval=t, rtol=1e-8, atol=1e-14)
    Y = s.y.reshape(16, 5, -1).transpose(2, 0, 1)
    ox = np.trapezoid(p["k450"] * Y[:, :, 0], s.t, axis=0) / p0
    gl = np.trapezoid(p["kG"] * Y[:, :, 0], s.t, axis=0) / p0
    return ox, gl, Y[:, :, 4].max(axis=0) / p0
schemes = [("Assumed (all)", "assumed", OR, "--"), ("Assumed, measured $k_G$ only", zm.perturb_assumed("kG"), GR, "-."),
           ("Measured (all)", "measured", BL, "-")]
fig, axes = plt.subplots(2, 3, figsize=(10.5, 6), sharex=True)
for r, dose in enumerate((4.0, 16.0)):
    res = [(lab, c, ls, sim(s, dose)) for lab, s, c, ls in schemes]
    for j, (title, k, scale) in enumerate([("Glucuronidated fraction of dose", 1, 100),
                                           ("Oxidised (NAPQI-forming) fraction of dose", 0, 100),
                                           ("Peak adduct fraction of dose", 2, 100)]):
        ax = axes[r, j]
        for lab, c, ls, v in res: ax.plot(x, v[k] * scale, color=c, ls=ls, lw=2.2, label=lab)
        ax.set_title(f"{title} (%), {dose:g} g", fontsize=10); ax.grid(alpha=0.3)
        if r == 1: ax.set_xlabel("hepatocyte (1 = periportal, 16 = pericentral)", fontsize=9.5)
axes[0, 0].legend(frameon=False, fontsize=9)
fig.tight_layout()
for ext in ("pdf", "svg"): fig.savefig(OUT / f"figS1_kG_mechanism.{ext}")
fig.savefig(OUT / "figS1_kG_mechanism.png", dpi=300)
print("done", sorted(p.name for p in OUT.iterdir()))
