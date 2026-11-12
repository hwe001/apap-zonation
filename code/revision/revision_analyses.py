"""Revision analyses for COMTOX-D-26-00340 (uses the repo code unmodified).

B4  hepatic penetration f sensitivity
B6  k_clear sensitivity
B8  range of isoform-weighted CYP fold across the 300 MC draws (same seed logic)
B11 mass balance by route + profiles, assumed vs measured kG (and full schemes)
"""
import sys, json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

REPO = Path(sys.argv[1]) / "code"
OUT = Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REPO))
import intracellular_apap_model as im
import zonated_apap_model as zm
import uncertainty as un

def peak_pericentral(dose, scheme):
    _, Y, _ = zm.run(scheme, dose, t_end_days=3.0, n_points=400, rtol=1e-5, atol=1e-8)
    return Y[:, :, 4].max(axis=0)[-1]

def change(dose, scheme_meas="measured"):
    b = peak_pericentral(dose, "assumed"); m = peak_pericentral(dose, scheme_meas)
    return (m / b - 1) * 100, b, m

res = {}

# ---- B4: f sensitivity ----
f0 = im.DOSAGE_PENETRATION
res["B4_f"] = {}
for f in (0.6, 0.8, 1.0):
    im.DOSAGE_PENETRATION = f
    res["B4_f"][f] = {d: change(d)[0] for d in (4.0, 8.0, 16.0)}
im.DOSAGE_PENETRATION = f0

# ---- B6: k_clear sensitivity ----
k0 = zm.K_CLEAR
res["B6_kclear"] = {}
for k in (0.2, 0.42, 1.0, 5.0):
    zm.K_CLEAR = k
    res["B6_kclear"][k] = {d: change(d)[0] for d in (4.0, 8.0, 16.0)}
zm.K_CLEAR = k0

# ---- B8: weighted CYP fold range over 300 draws (matches uncertainty.py seed) ----
rng = np.random.default_rng(0)
folds = np.array([un.sample_cyp_fold(rng) for _ in range(300)])
res["B8_cyp_fold"] = {"min": folds.min(), "p5": np.percentile(folds, 5),
                      "median": np.median(folds), "p95": np.percentile(folds, 95), "max": folds.max()}

# ---- B11: mass balance ----
def routes(scheme, dose):
    p = zm.build_parameters(scheme)
    p0 = im.dose_to_p0(dose); s0 = im.B_S / im.D_S; g0 = p["bG"] / im.D_G
    n = zm.N_HEPATOCYTES
    y0 = np.stack([np.full(n, p0), np.full(n, s0), np.zeros(n), g0, np.zeros(n)], 1).ravel()
    t = np.concatenate([[0], np.geomspace(1e-5, 3.0, 4000)])
    sol = solve_ivp(zm.make_rhs(p), (0, 3.0), y0, method="Radau", t_eval=t, rtol=1e-8, atol=1e-14)
    Y = sol.y.reshape(n, 5, -1).transpose(2, 0, 1)
    P, S = Y[:, :, 0], Y[:, :, 1]
    flux = {"glucuronidation": p["kG"] * P, "sulfation": p["kS"] * S * P, "oxidation": p["k450"] * P}
    tot = {k: np.trapezoid(v, sol.t, axis=0) for k, v in flux.items()}  # per cell, nmol
    return tot, p0, Y[:, :, 4].max(axis=0), Y[:, :, 0]

mb = {}
for dose in (4.0, 16.0):
    mb[dose] = {}
    for name, sch in {"assumed": "assumed", "measured": "measured",
                      "assumed_but_measured_kG": zm.perturb_assumed("kG")}.items():
        tot, p0, Cpk, P = routes(sch, dose)
        d = {"p0": p0, "peak_C": Cpk.tolist()}
        for r, v in tot.items():
            d[r] = {"pp_cell1": v[0] / p0, "pc_cell16": v[-1] / p0, "mean": v.mean() / p0}
        mb[dose][name] = d
res["B11"] = mb

def conv(o):
    if isinstance(o, (np.floating, np.integer)): return float(o)
    if isinstance(o, dict): return {str(k): conv(v) for k, v in o.items()}
    if isinstance(o, list): return [conv(x) for x in o]
    return o
json.dump(conv(res), open(OUT / "revision_results.json", "w"), indent=1)

print("B4 f sensitivity (% change measured vs assumed):")
for f, v in res["B4_f"].items(): print(" f=%.1f" % f, {k: round(x, 1) for k, x in v.items()})
print("B6 k_clear sensitivity:")
for k, v in res["B6_kclear"].items(): print(" k=%.2f" % k, {d: round(x, 1) for d, x in v.items()})
print("B8 weighted CYP fold:", {k: round(float(v), 2) for k, v in res["B8_cyp_fold"].items()})
print("B11 route fractions (fraction of initial dose per cell; PP=cell1, PC=cell16):")
for dose, sc in mb.items():
    for name, d in sc.items():
        print(f" {dose:g} g {name:25s}",
              " | ".join(f"{r[:4]} PP {d[r]['pp_cell1']:.3f} PC {d[r]['pc_cell16']:.3f}" for r in ("glucuronidation", "sulfation", "oxidation")),
              f"| peakC_PC {d['peak_C'][-1]:.3e}")
