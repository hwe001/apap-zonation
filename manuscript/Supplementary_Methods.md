# Supplementary Methods — fold derivation from the DVP bin profiles

## S1. Source and scaling

All protein-level data are from Weiss et al. (Nat Metab 2026), Supplementary
Table 3 (healthy cohort, N = 14, proteins with >=70% data completeness). The table
reports, per protein: (i) mean expression intensity in 20 equidistant bins along
the porto-central axis, normalised so that the 20 bin values sum to 100 (bin 0 =
central-vein end, bin 19 = portal-vein end); and (ii) a "zonation coefficient",
defined as the fixed-effect slope of a per-protein linear mixed model fitted to
robust-scaled (median-centred, IQR-divided) intensities. The coefficient is a
regression slope on scaled data and is NOT a portal/central abundance ratio; it is
reported below for comparison only and is not used in the fold derivation.

## S2. Fold calculation

For each protein:
1. C = mean of bins 0-4 (central-most fifth of the sinusoid);
2. P = mean of bins 15-19 (portal-most fifth);
3. fold = C / P (a value > 1 means central-enriched).

The fold is applied in the model as a pericentral/periportal rate multiplier with
an arithmetic-mean-preserving linear ramp (Methods 2.2), so total enzyme is
conserved. Using five bins per end rather than the single terminal bin reduces
sensitivity to noise in any one bin; using the extreme bins alone changes the
folds by <30% (single-bin ratios are listed below) and does not change any
conclusion.

## S3. Per-protein bin values (mean intensity per bin, % of total)

Bins run central (0) to portal (19). The proteins mapped to model parameters in
Table 2 of the main text are CYP2E1/CYP3A4/CYP1A2 (k450), UGT2B7/UGT1A1 (kG),
SULT1A1 (kS), GSTA2 (kGSH), and GCLC/GCLM/GSS (bG); the remaining rows are
landmark and validation proteins shown for context.

| protein | b0 | b1 | b2 | b3 | b4 | b5 | b6 | b7 | b8 | b9 | b10 | b11 | b12 | b13 | b14 | b15 | b16 | b17 | b18 | b19 | C (0-4) | P (15-19) | fold C/P | single-bin b0/b19 | coefficient |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| UROC1 | 2.46 | 3.23 | 3.26 | 3.32 | 3.87 | 4.01 | 3.41 | 4.76 | 5.47 | 5.08 | 5.52 | 5.96 | 5.76 | 6.41 | 5.91 | 6.68 | 5.61 | 7.04 | 6.39 | 5.86 | 3.23 | 6.32 | 0.51 | 0.42 | 2.09 |
| AMDHD1 | 3.08 | 2.93 | 2.95 | 2.20 | 3.03 | 4.21 | 2.89 | 4.03 | 4.68 | 5.23 | 5.80 | 5.11 | 5.72 | 6.31 | 6.08 | 7.36 | 5.72 | 7.53 | 6.81 | 8.34 | 2.84 | 7.15 | 0.40 | 0.37 | 1.98 |
| ASS1 | 2.72 | 3.10 | 2.99 | 2.91 | 3.32 | 3.97 | 3.44 | 5.03 | 4.37 | 5.48 | 5.15 | 6.12 | 5.76 | 5.63 | 6.21 | 7.63 | 6.23 | 7.00 | 6.38 | 6.56 | 3.01 | 6.76 | 0.44 | 0.42 | 1.89 |
| CPS1 | 3.54 | 4.17 | 3.99 | 4.26 | 4.54 | 4.62 | 3.98 | 5.51 | 4.94 | 4.87 | 5.23 | 5.39 | 5.36 | 5.57 | 5.42 | 5.89 | 5.61 | 5.63 | 5.98 | 5.49 | 4.10 | 5.72 | 0.72 | 0.64 | 1.79 |
| SULT1A1 | 4.28 | 5.15 | 5.29 | 4.55 | 5.12 | 4.44 | 5.41 | 4.91 | 5.51 | 5.64 | 5.44 | 5.47 | 5.14 | 4.41 | 4.67 | 5.97 | 4.34 | 4.64 | 4.87 | 4.79 | 4.88 | 4.92 | 0.99 | 0.89 | 0.25 |
| GCLM | 5.44 | 5.02 | 5.02 | 4.36 | 5.19 | 4.85 | 4.92 | 4.71 | 4.72 | 5.36 | 5.13 | 4.68 | 5.18 | 4.78 | 4.79 | 4.95 | 5.30 | 5.02 | 5.35 | 5.24 | 5.00 | 5.17 | 0.97 | 1.04 | 0.16 |
| GCLC | 5.19 | 4.82 | 4.61 | 5.00 | 5.22 | 5.19 | 5.15 | 4.99 | 4.92 | 5.38 | 5.41 | 4.51 | 5.10 | 5.03 | 5.08 | 5.36 | 4.81 | 4.76 | 4.87 | 4.60 | 4.97 | 4.88 | 1.02 | 1.13 | -0.05 |
| GSR | 5.38 | 5.07 | 4.59 | 5.33 | 4.51 | 4.90 | 5.57 | 4.72 | 5.22 | 4.59 | 5.11 | 4.79 | 5.08 | 5.07 | 5.19 | 4.37 | 4.91 | 4.75 | 5.54 | 5.31 | 4.98 | 4.98 | 1.00 | 1.01 | -0.08 |
| GSS | 5.59 | 5.69 | 5.20 | 4.83 | 4.52 | 5.09 | 5.00 | 5.19 | 4.81 | 4.90 | 5.14 | 4.42 | 5.31 | 4.13 | 4.92 | 4.45 | 4.62 | 5.89 | 5.03 | 5.28 | 5.16 | 5.06 | 1.02 | 1.06 | -0.08 |
| GSTT1 | 5.38 | 5.60 | 4.72 | 4.62 | 5.36 | 4.61 | 5.64 | 4.51 | 5.03 | 4.88 | 6.06 | 5.46 | 4.84 | 5.38 | 4.23 | 5.55 | 4.25 | 4.85 | 4.27 | 4.78 | 5.13 | 4.74 | 1.08 | 1.12 | -0.22 |
| GPX1 | 5.09 | 5.88 | 5.85 | 5.21 | 4.24 | 4.90 | 4.88 | 4.89 | 5.04 | 4.89 | 4.50 | 5.07 | 5.22 | 5.61 | 4.82 | 4.94 | 4.73 | 4.42 | 5.35 | 4.47 | 5.25 | 4.78 | 1.10 | 1.14 | -0.51 |
| GSTA2 | 5.07 | 4.97 | 5.18 | 5.50 | 5.10 | 5.02 | 6.65 | 4.84 | 5.22 | 5.09 | 4.80 | 5.03 | 5.12 | 4.82 | 5.35 | 4.17 | 4.81 | 4.53 | 4.50 | 4.24 | 5.16 | 4.45 | 1.16 | 1.20 | -0.60 |
| UGT1A3 | 5.89 | 6.83 | 6.43 | 7.44 | 4.78 | 6.90 | 8.41 | 5.72 | 3.45 | 6.49 | 4.24 | 3.67 | 5.58 | 1.75 | 5.77 | 3.71 | 4.14 | 2.40 | 3.14 | 3.25 | 6.27 | 3.33 | 1.88 | 1.81 | -0.85 |
| UGT1A4 | 6.64 | 6.23 | 6.27 | 6.63 | 4.96 | 6.30 | 6.48 | 5.37 | 3.95 | 5.52 | 5.41 | 4.09 | 5.36 | 3.27 | 4.95 | 3.91 | 4.93 | 2.69 | 3.33 | 3.72 | 6.14 | 3.72 | 1.65 | 1.78 | -0.95 |
| UGT1A9 | 6.56 | 5.67 | 6.20 | 6.30 | 5.40 | 5.48 | 4.99 | 4.38 | 4.55 | 4.79 | 4.45 | 4.65 | 5.25 | 4.63 | 4.49 | 4.55 | 4.00 | 3.89 | 5.03 | 4.74 | 6.03 | 4.44 | 1.36 | 1.38 | -1.06 |
| UGT2B7 | 6.08 | 5.34 | 5.73 | 6.43 | 5.26 | 5.56 | 6.04 | 4.75 | 4.86 | 5.79 | 5.20 | 4.62 | 5.09 | 4.50 | 4.79 | 3.90 | 4.56 | 3.77 | 4.10 | 3.63 | 5.77 | 3.99 | 1.44 | 1.67 | -1.09 |
| SULT2A1 | 5.86 | 6.57 | 5.68 | 6.45 | 5.69 | 5.69 | 6.77 | 5.39 | 4.70 | 5.50 | 6.01 | 4.86 | 5.14 | 4.02 | 4.87 | 3.97 | 4.32 | 2.75 | 3.27 | 2.51 | 6.05 | 3.36 | 1.80 | 2.33 | -1.20 |
| UGT1A1 | 7.88 | 5.85 | 6.96 | 7.51 | 6.09 | 5.46 | 6.39 | 4.47 | 5.32 | 5.40 | 4.86 | 4.47 | 5.10 | 3.39 | 4.73 | 2.99 | 4.82 | 1.84 | 3.38 | 3.09 | 6.86 | 3.23 | 2.13 | 2.55 | -1.24 |
| CYP3A4 | 8.86 | 7.83 | 7.12 | 8.10 | 7.06 | 6.56 | 9.02 | 5.05 | 3.21 | 7.07 | 5.51 | 3.37 | 4.83 | 1.97 | 3.38 | 2.73 | 4.16 | 0.88 | 1.43 | 1.86 | 7.79 | 2.21 | 3.52 | 4.76 | -1.60 |
| UGT1A6 | 7.52 | 7.28 | 7.53 | 6.19 | 5.93 | 6.39 | 5.46 | 4.53 | 4.66 | 4.65 | 4.61 | 4.16 | 4.61 | 3.79 | 3.78 | 3.98 | 4.09 | 2.79 | 4.07 | 3.98 | 6.89 | 3.78 | 1.82 | 1.89 | -1.68 |
| CYP2E1 | 6.07 | 5.78 | 6.48 | 6.74 | 5.95 | 5.83 | 6.18 | 5.98 | 5.53 | 5.68 | 4.92 | 5.03 | 4.99 | 4.52 | 4.05 | 3.80 | 4.17 | 2.69 | 3.07 | 2.55 | 6.20 | 3.25 | 1.91 | 2.38 | -1.82 |
| CYP1A2 | 12.80 | 11.90 | 8.69 | 9.60 | 7.54 | 5.22 | 8.64 | 4.59 | 3.91 | 5.36 | 4.93 | 2.60 | 3.24 | 1.59 | 2.03 | 1.46 | 1.67 | 1.02 | 1.54 | 1.66 | 10.11 | 1.47 | 6.87 | 7.70 | -2.20 |

## S4. Mapping to model parameters and uncertainty inputs

Table 2 of the main text maps proteins to the five model parameters. Where
multiple isoforms contribute, the model uses one representative isoform (stated
in the main text) and the isoform span defines the uncertainty range. The Monte
Carlo samples each fold log-uniformly over its range (400 draws):

| parameter | isoform used | point fold | MC range | basis of range |
|---|---|---|---|---|
| k450 (CYP450) | CYP2E1/1A2/3A4 weighted 0.60/0.25/0.15 | 3.39 | sampled via weight vector (2E1 0.50-0.70, 1A2 0.15-0.30, 3A4 remainder) | weights = approximate relative contribution to human APAP oxidation (Laine et al. 2009); raw isoform folds 1.91/6.87/3.52 |
| kG (UGT) | UGT2B7 | 1.44 | 1.4 - 2.1 | isoform span UGT2B7..UGT1A1 |
| kS (SULT) | SULT1A1 | 0.99 | 0.9 - 1.1 | not significant (q=0.12); +/-10% around uniform (stated assumption) |
| kGSH (GST) | GSTA2 | 1.16 | 1.16 - 3.4 | protein fold .. mRNA ratio (Yakubovsky 2026) |
| bG (GSH synthesis) | GCLC/GCLM/GSS | ~1.0 | 0.9 - 1.1 | coefficients ~0; +/-10% around uniform (stated assumption) |

Donor-level variation is not propagated at this stage: the bin profiles are
cohort-level means across the 14 healthy donors, and the published supplement
does not report per-donor bin values. The simulation intervals in Fig. 4
therefore represent scenario uncertainty over fold ranges, not sampling
uncertainty estimated from donors.

## S5. Validation landmarks

The same calculation applied to established landmarks reproduces the expected
directions, supporting the pipeline: ASS1, CPS1, UROC1 and AMDHD1 (periportal
markers of the urea cycle and histidine catabolism) yield folds < 1
(portal-enriched), while CYP1A2, CYP3A4 and CYP2E1 yield folds > 1
(central-enriched), matching their established zonation.
- ASS1: fold = 0.44 (coefficient 1.89)
- CPS1: fold = 0.72 (coefficient 1.79)
- UROC1: fold = 0.51 (coefficient 2.09)
- AMDHD1: fold = 0.40 (coefficient 1.98)
- CYP1A2: fold = 6.87 (coefficient -2.20)
- CYP3A4: fold = 3.52 (coefficient -1.60)
