"""Quick checks: (1) time of peak pericentral adduct, (2) which MC fold explains 8 g -31% vs -37%,
(3) Fig. 4 re-plot from the reproduced Monte Carlo output."""
import sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(sys.argv[1]); OUT = Path(sys.argv[2])
sys.path.insert(0, str(REPO / "code"))
import zonated_apap_model as zm

def pk(scheme, dose):
    _, Y, _ = zm.run(scheme, dose, t_end_days=3.0, n_points=400, rtol=1e-5, atol=1e-8)
    return Y[:, :, 4].max(axis=0)[-1]

# (1) time of peak
print("(1) time of peak pericentral adduct (d):")
for dose in (4.0, 8.0, 16.0):
    for sch in ("assumed", "measured"):
        t, Y, _ = zm.run(sch, dose, t_end_days=3.0, n_points=3000, rtol=1e-8, atol=1e-14)
        print(f"   {dose:g} g {sch:9s}: t_peak = {t[np.argmax(Y[:, -1, 4])]:.3f} d")

# (2) 8 g decomposition: point scheme, then each sampled fold moved to the MC median / ends
pt = {"k450": 3.39, "kG": 1.44, "kS": 1.0, "kGSH": 1.16, "bG": 1.0}
mc_mid = {"k450": 3.29, "kG": float(np.sqrt(1.4 * 2.1)), "kS": 1.0, "kGSH": float(np.sqrt(1.16 * 3.4)), "bG": 1.0}
def scheme(f):
    return {p: ((1.0, 1.0) if f[p] == 1.0 and p in ("kS", "bG") else zm.am_fold(f[p], True)) for p in f}
print("(2) 8 g change vs assumed; point folds, then one fold at its MC log-midpoint:")
for dose in (8.0,):
    base = pk("assumed", dose)
    print(f"   point folds           : {(pk(scheme(pt), dose)/base-1)*100:+.1f}%")
    for p in ("k450", "kG", "kGSH"):
        f = dict(pt); f[p] = mc_mid[p]
        print(f"   {p:4s} -> {mc_mid[p]:.2f}          : {(pk(scheme(f), dose)/base-1)*100:+.1f}%")
    print(f"   all three at MC mid   : {(pk(scheme(mc_mid), dose)/base-1)*100:+.1f}%")

# (3) Fig. 4
d = np.load(REPO / "results" / "uncertainty_changes.npz")
plt.rcParams.update({"font.size": 11, "pdf.fonttype": 42})
fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
for ax, dose, key in zip(axes, (4, 8, 16), ("changes_4g", "changes_8g", "changes_16g")):
    c = d[key]; lo, med, hi = np.percentile(c, [5, 50, 95])
    ax.hist(c, bins=24, color="#2a78d6", alpha=0.8, edgecolor="white")
    ax.axvline(0, color="0.3", lw=1); ax.axvline(med, color="#eb6834", lw=2, ls="--")
    ax.set_title(f"{dose} g single dose"); ax.set_xlabel("change in peak pericentral adducts (%)")
    ax.text(0.97, 0.95, f"median {med:+.1f}%\n90% interval [{lo:+.1f}, {hi:+.1f}]\n{(c<0).mean()*100:.0f}% protective, {(c>0).mean()*100:.0f}% aggravating",
            transform=ax.transAxes, ha="right", va="top", fontsize=9.5)
axes[0].set_ylabel("Monte Carlo draws (of 300)")
fig.tight_layout()
for ext in ("pdf", "svg"): fig.savefig(OUT / f"fig4_uncertainty_revised.{ext}")
fig.savefig(OUT / "fig4_uncertainty_revised.png", dpi=600)
print("(3) wrote fig4_uncertainty_revised")
