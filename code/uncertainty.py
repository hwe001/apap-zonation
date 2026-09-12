#!/usr/bin/env python3
"""Uncertainty propagation for the zonation recalibration (JTB revision).

Monte Carlo over the defensible fold ranges for each enzyme — the isoform span
plus the protein/mRNA discordance — to test whether the claim "assumed
gradients overestimate pericentral adduct burden" survives when the measured
folds are uncertain. The assumed (literature) gradients are held fixed as the
null hypothesis being tested; only the measured folds are sampled.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from zonated_apap_model import run, am_fold  # noqa: E402

# Fold ranges (central/portal), defensible spans:
RANGES = {
    "k450": (1.9, 6.9),    # CYP2E1 .. CYP1A2 protein folds (CYP3A4 3.5x inside)
    "kG":   (1.4, 2.1),    # UGT2B7 .. UGT1A1 protein folds
    "kS":   (0.9, 1.1),    # SULT1A1 ~uniform, +/-10% (uncertainty on "unzonated")
    "kGSH": (1.16, 3.4),   # GSTA2 protein (1.16x) .. mRNA (3.4x)
    "bG":   (0.9, 1.1),    # GCLC/GCLM/GSS ~uniform, +/-10%
}


def sample_measured(rng):
    """Sample a 'measured' scheme: each fold drawn log-uniform over its range."""
    s = {}
    for p, (lo, hi) in RANGES.items():
        F = np.exp(rng.uniform(np.log(lo), np.log(hi)))
        s[p] = am_fold(F, central=True)  # F~1 for kS/bG -> near-uniform either way
    return s


def pericentral_peak(dose, scheme):
    _, Y, _ = run(scheme, dose, t_end_days=3.0, n_points=400, rtol=1e-5, atol=1e-8)
    return Y[:, :, 4].max(axis=0)[-1]  # peak adduct, pericentral cell


def main():
    rng = np.random.default_rng(0)
    N = 400
    saved = {}

    for dose in (4.0, 16.0):
        base = pericentral_peak(dose, "assumed")
        changes = []
        for _ in range(N):
            s = sample_measured(rng)
            c = pericentral_peak(dose, s)
            changes.append((c / base - 1) * 100)
        changes = np.array(changes)
        saved[dose] = changes
        lo, med, hi = np.percentile(changes, [5, 50, 95])
        print(f"{dose:.0f} g  ({N} draws):")
        print(f"  change in pericentral adducts  median {med:+.1f}%  "
              f"90% CI [{lo:+.1f}, {hi:+.1f}]")
        print(f"  fraction protective (reduction): {np.mean(changes < 0) * 100:.0f}%")
        print(f"  fraction aggravating (increase): {np.mean(changes > 0) * 100:.0f}%")
        print()

    out = Path(__file__).resolve().parent.parent / "results" / "uncertainty_changes.npz"
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(out, changes_4g=saved[4.0], changes_16g=saved[16.0])
    print(f"saved {out}")


if __name__ == "__main__":
    main()
