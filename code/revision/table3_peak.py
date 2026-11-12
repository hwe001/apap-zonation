"""Recompute Table 3 / point estimates with (a) peak-over-time metric (as in uncertainty.py) and
(b) end-of-simulation (t=5 d) value (as in zonated_apap_model.decompose), to expose the metric mismatch."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(sys.argv[1]) / "code"))
import zonated_apap_model as zm
import intracellular_apap_model as im

def metrics(scheme, dose):
    t, Y, _ = zm.run(scheme, dose, t_end_days=5.0, n_points=1200, rtol=1e-8, atol=1e-14)
    C = Y[:, -1, 4]
    return C.max(), C[-1]   # peak (pericentral cell), value at 5 d

for dose in (4.0, 8.0, 16.0):
    bpk, bend = metrics("assumed", dose); p0 = im.dose_to_p0(dose)
    print(f"{dose:g} g assumed: peak {bpk:.3e} ({bpk/p0*100:.2f}% dose) | day-5 {bend:.3e} ({bend/p0*100:.2f}% dose)")
    for name, sch in [("ALL measured", "measured")] + [(p, zm.perturb_assumed(p)) for p in ("kG", "bG", "k450", "kS", "kGSH")]:
        pk, en = metrics(sch, dose)
        print(f"   {name:13s} peak {((pk/bpk)-1)*100:+6.1f}%   day-5 {((en/bend)-1)*100:+6.1f}%")
