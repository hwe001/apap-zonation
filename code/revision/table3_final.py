"""Corrected Table 3, point estimates and 27-variant null sensitivity, all on the TRUE PEAK
pericentral adduct (max over time of cell 16; t_end 5 d, tight tolerances)."""
import sys, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(sys.argv[1]) / "code"))
import zonated_apap_model as zm
import intracellular_apap_model as im

def pk(scheme, dose):
    _, Y, _ = zm.run(scheme, dose, t_end_days=5.0, n_points=1200, rtol=1e-8, atol=1e-14)
    return Y[:, -1, 4].max()

out = {"point": {}, "table3": {}, "null27": {}}
for dose in (4.0, 8.0, 16.0):
    p0 = im.dose_to_p0(dose); a = pk("assumed", dose); m = pk("measured", dose)
    out["point"][dose] = dict(assumed=a, measured=m, assumed_pct=a/p0*100, measured_pct=m/p0*100, change=(m/a-1)*100)
    out["table3"][dose] = {p: (pk(zm.perturb_assumed(p), dose)/a-1)*100 for p in ("kG", "bG", "k450", "kS", "kGSH")}
m4 = pk("measured", 4.0); d = []
for fc in (3., 5., 7.):
    for fu in (6., 9., 12.):
        for fp in (6., 9., 12.):
            s = {"k450": zm.am_fold(fc, True), "kG": zm.am_fold(fu, True),
                 "kS": zm.am_fold(fp, False), "kGSH": zm.am_fold(fp, False), "bG": zm.am_fold(fp, False)}
            d.append((m4/pk(s, 4.0)-1)*100)
d = np.array(d)
out["null27"] = dict(median=float(np.median(d)), min=float(d.min()), max=float(d.max()), frac_protective=float((d < 0).mean()))
json.dump(out, open(Path(sys.argv[2]) / "table3_peak_final.json", "w"), indent=1, default=float)
for dose in (4.0, 8.0, 16.0):
    pt = out["point"][dose]
    print(f"{dose:g} g: assumed {pt['assumed']:.2e} nmol ({pt['assumed_pct']:.2f}% dose) -> measured {pt['measured']:.2e} ({pt['measured_pct']:.2f}%)  change {pt['change']:+.1f}%")
    print("   OAT:", {k: round(v, 1) for k, v in out["table3"][dose].items()})
print("null27 (4 g, peak):", out["null27"])
