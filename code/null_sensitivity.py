#!/usr/bin/env python3
"""Sensitivity of the stability result to the exact values of the ASSUMED (null)
gradients.

The prior models impose qualitative "strong" gradients rather than precise
folds; the 5x/9x values used in the main text are stylised. This script varies
the assumed folds over plausible stylised ranges (CYP 3-7x, UGT and the
periportal gradients 6-12x), holds the measured gradients fixed, and reports
the resulting change in pericentral adduct burden at 4 g for every combination.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from zonated_apap_model import SCHEMES, run, am_fold, summarize  # noqa: E402

DOSE = 4.0
CYP_FOLDS = (3.0, 5.0, 7.0)
UGT_FOLDS = (6.0, 9.0, 12.0)
PORTAL_FOLDS = (6.0, 9.0, 12.0)  # kS, kGSH, bG assumed periportal folds


def assumed_scheme(F_cyp, F_ugt, F_portal):
    return {
        "k450": am_fold(F_cyp, True),
        "kG": am_fold(F_ugt, True),
        "kS": am_fold(F_portal, False),
        "kGSH": am_fold(F_portal, False),
        "bG": am_fold(F_portal, False),
    }


def main():
    _, Y_m, _ = run("measured", DOSE)
    meas = summarize(None, Y_m)[2]["central_C"]
    print(f"measured pericentral C (4 g) = {meas:.4g} nmol\n")

    print(f"{'CYP':>4s} {'UGT':>4s} {'portal':>7s} | {'assumed C':>10s} {'measured/assumed-1':>19s}")
    print("-" * 52)
    deltas = []
    for F_cyp in CYP_FOLDS:
        for F_ugt in UGT_FOLDS:
            for F_portal in PORTAL_FOLDS:
                _, Y_a, _ = run(assumed_scheme(F_cyp, F_ugt, F_portal), DOSE)
                a = summarize(None, Y_a)[2]["central_C"]
                d = (meas / a - 1) * 100
                deltas.append(d)
                print(f"{F_cyp:4.0f} {F_ugt:4.0f} {F_portal:7.0f} | {a:10.4g} {d:+18.1f}%")
    deltas = np.array(deltas)
    print(f"\nall {len(deltas)} assumed-gradient variants: "
          f"median {np.median(deltas):+.1f}%, range [{deltas.min():+.1f}, {deltas.max():+.1f}]%")
    print(f"protective in {np.mean(deltas < 0) * 100:.0f}% of variants")


if __name__ == "__main__":
    main()
