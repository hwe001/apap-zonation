# Changelog

## v1.2.0 — revision (COMTOX-D-26-00340)

- `zonated_apap_model.py`: `decompose()` now reports the **peak** pericentral adduct (it previously used the t = 5 d value, which did not match the Monte Carlo or Fig. 3). One-at-a-time effects at the peak: 4 g kG +41.6%, bG −83.6%, k450 −8.6%, kS −8.4%, kGSH −2.2%; all measured −90.0%.
- Note: the summary lines printed by `main()` (`central_C`, "kGSH flip alone", "full recalibration") still use the t = 5 d value. Use the one-at-a-time block or `code/revision/table3_final.py` for peak-based numbers.
- `zonated_apap_model.py`: docstring fold values corrected to match the code (k450 3.39×, kG 1.44×, kGSH 1.16×).
- `requirements.txt`: versions pinned (numpy 2.4.6, scipy 1.17.1, matplotlib 3.11.2; Python 3.11.9).
- New `code/revision/`: sensitivity to hepatic penetration f and adduct clearance k_clear, weighted-CYP-fold range, glucuronidation mass balance, corrected Table 3, peak-time and 8 g checks, revised Fig. 2/4 and Fig. S1. Run order: `uncertainty.py` first (writes `results/uncertainty_changes.npz`), then the scripts in `code/revision/` with `<repo root> <output dir>` arguments.
- New `figures/revision/`: revised figures and result files.
- Monte Carlo (seed 0, 300 draws) re-run: 4 g −93.5% (−95.6 to −88.9), 8 g −36.9% (−43.2 to −30.0), 16 g −6.9% (−12.8 to +0.3), 94% protective at 16 g. Unchanged from the first submission.

## v1.1.1 — first submission
