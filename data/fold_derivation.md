# Fold derivation — per-protein inputs

## Method

From Weiss et al. (Nat Metab 2026), Supplementary Table 3, which reports each
protein's mean intensity in 20 equidistant bins along the porto-central axis,
normalised so the 20 bins sum to 100% (bin 0 = central vein, bin 19 = portal vein).
For each protein we:
1. average the five central-most bins (bins 0–4);
2. average the five portal-most bins (bins 15–19);
3. take the central/portal ratio of those two averages.

The published "zonation coefficient" (the robust-scaled LMM slope) is shown for
comparison only — it is not a fold change and is not used.

## Per-protein inputs (mean intensity per bin, % of total)

| gene | bins 0–4 (central) | bins 15–19 (portal) | central/portal fold | coefficient |
|---|---|---|---|---|
| ASS1 | 3.01 | 6.76 | 0.44 | 1.89 |
| CPS1 | 4.10 | 5.72 | 0.72 | 1.79 |
| SULT1A1 | 4.88 | 4.92 | 0.99 | 0.25 |
| GCLM | 5.00 | 5.17 | 0.97 | 0.16 |
| GCLC | 4.97 | 4.88 | 1.02 | -0.05 |
| GSR | 4.98 | 4.98 | 1.00 | -0.08 |
| GSS | 5.16 | 5.06 | 1.02 | -0.08 |
| GPX1 | 5.25 | 4.78 | 1.10 | -0.51 |
| GSTA2 | 5.16 | 4.45 | 1.16 | -0.60 |
| UGT1A9 | 6.03 | 4.44 | 1.36 | -1.06 |
| UGT2B7 | 5.77 | 3.99 | 1.44 | -1.09 |
| SULT2A1 | 6.05 | 3.36 | 1.80 | -1.20 |
| UGT1A1 | 6.86 | 3.23 | 2.13 | -1.24 |
| CYP3A4 | 7.79 | 2.21 | 3.52 | -1.60 |
| UGT1A6 | 6.89 | 3.78 | 1.82 | -1.68 |
| CYP2E1 | 6.20 | 3.25 | 1.91 | -1.82 |
| CYP1A2 | 10.11 | 1.47 | 6.87 | -2.20 |
