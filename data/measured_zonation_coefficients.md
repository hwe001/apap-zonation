# Measured zonal zonation coefficients — extracted from the two 2026 human datasets

**Date compiled:** 2026-09-13
**Purpose:** Real numbers to plug into the 2018 model's zonal gradients (kS/kG/k450/kGSH/bG), replacing the qualitative Gebhardt (1992) gradient. Extracted from the papers' supplementary tables (downloaded and parsed directly, not the main-text summaries in `spatial_liver_zonation_datasets.md`).

## Two datasets, two levels

| Source | Species | Level | What the number is | File parsed |
|---|---|---|---|---|
| **DVP** — *Nat Metab* 2026 ([PMC13031132](https://pmc.ncbi.nlm.nih.gov/articles/PMC13031132/), DOI 10.1038/s42255-026-01459-2) | Human, N=14 healthy | **Protein** | "zonation coefficient" (regression), Supplementary Table 3 | `MOESM3_ESM.xlsx`, sheet "Table 3" |
| **Yakubovsky** — *Nature* 2026 ([PMC13216088](https://pmc.ncbi.nlm.nih.gov/articles/PMC13216088/), DOI 10.1038/s41586-026-10377-y) | Human (+4 mammals) | **mRNA** (Visium) | log2(portal/central) ratio, Fig. 3a source data | `MOESM5_ESM.xlsx`, sheet "Fig3a" |

**Direction convention (DVP coefficient):** sign carries direction — **negative = pericentral (zone 3), positive = periportal (zone 1)**. Verified against canonical landmarks in the same table: ASS1 (+1.89) and CPS1 (+1.79) are periportal; CYP3A4/CYP2E1/CYP1A2 are pericentral. (Bins run pericentral→periportal; bin 0 = central.)

For the model, the **protein-level (DVP) numbers are the more directly relevant** — the model parameters are reaction rates, i.e. enzyme activity, not transcript. mRNA (Yakubovsky) is reported as a cross-check; the two agree on every direction.

---

## CYP450 → k450 (model H4)

| Gene | DVP coeff | Yakubovsky log2(portal/central) human / mouse | Direction |
|---|---|---|---|
| CYP3A4 | **−1.60** | — (no mouse orthologue) | pericentral |
| CYP2E1 | **−1.82** | −2.42 / −3.11 | pericentral |
| CYP1A2 | **−2.20** | −3.50 / −2.72 | pericentral |

All three pericentral, all strongly zonated (q ≤ 1e-33). Confirms H4's direction and supplies quantitative magnitudes. CYP1A2 is the most strongly zonated of the three.

## Glucuronidation → kG (model H3)

| Gene | DVP coeff | Direction |
|---|---|---|
| UGT2B7 | **−1.09** | pericentral |
| UGT1A6 | −1.68 | pericentral |
| UGT1A1 | −1.24 | pericentral |
| UGT1A9 | −1.06 | pericentral |
| UGT1A4 | −0.95 | pericentral |
| UGT1A3 | −0.85 | pericentral |

Every detectable UGT is pericentral. Confirms H3's direction with protein data. Note UGT2B7 has no clean mouse orthologue, so it is absent from Yakubovsky's cross-species table — the DVP protein value is the one to use.

## Sulfation → kS (model H2)

| Gene | DVP coeff | zonated? | Direction |
|---|---|---|---|
| **SULT1A1** (the APAP sulfotransferase) | **+0.25** | **No** (q=0.12) | **~unzonated** |
| SULT2A1 (DHEA-sulfating, not APAP-relevant) | −1.20 | Yes | pericentral |

Yakubovsky mRNA agrees: SULT1A1 log2(portal/central) = −0.45 (human) vs +0.02 (mouse) — weak and inconsistent in sign. **The enzyme that actually sulfates APAP is not significantly zonated.** This undercuts H2's "sulfation higher periportal" assumption — and matches the equivocal literature (Anundi 1993 found no periportal dominance; MALDI-MSI found no zonation of the sulfate conjugate).

## Glutathione → kGSH (H6, conjugation) and bG/G0 (H5, synthesis)

| Gene | DVP coeff | zonated? | Yakubovsky log2 human / mouse | Direction |
|---|---|---|---|---|
| **GSTA2** (→ kGSH) | **−0.60** | Yes | **−1.75** / −0.17 | **pericentral** |
| GSTT1 | −0.22 | No | — | weak pericentral |
| GPX1 | −0.51 | Yes | — | pericentral |
| GCLC (→ bG) | −0.05 | No | −0.68 / −0.42 | ~unzonated |
| GCLM (→ bG) | +0.16 | No | −0.18 / −0.40 | ~unzonated |
| GSS (→ bG) | −0.08 | No | −0.15 / +0.47 | ~unzonated |
| GSR | −0.08 | No | 0.00 / −0.51 | ~unzonated |

(GSTA1, GSTM1, GSTP1 not detected in the DVP panel — below detection limit.)

Two independent findings:
1. **GST conjugation (kGSH) is pericentral** — DVP protein (−0.60), Yakubovsky mRNA (−1.75), and total GST *activity* (El Mouelhi & Kauffman 1986, human: 12.2 pericentral vs 8.3 periportal). This is the opposite of the model's H6 ("GSH binding higher periportal").
2. **GSH synthesis (GCLC/GCLM/GSS) is unzonated** — the model's H5 ("GSH production higher periportal") is not supported by the enzyme distribution (GSH *levels* may still be higher periportal for consumption/transport reasons, but the synthetic enzymes are uniform).

---

## Payoff: model parameter → measured → action

| Model param | Prior papers' assumption | Measured (human) | Action for recalibration |
|---|---|---|---|
| **k450** (CYP450) | pericentral, nonlinear/IHC | CYP3A4 −1.60, CYP2E1 −1.82, CYP1A2 −2.20 | **Confirm + make quantitative** (isoform-specific) |
| **kG** (UGT) | pericentral | UGT2B7 −1.09 | **Confirm + quantitative** |
| **kS** (SULT) | periportal | SULT1A1 +0.25, **unzonated** | **Remove gradient** (near-uniform) |
| **kGSH** (GST) | periportal | GSTA2 −0.60/−1.75, **pericentral** | **Sign flip** — reverse the gradient |
| **bG/G0** (GSH synth) | periportal | GCLC/GCLM/GSS ~0, **unzonated** | **Remove gradient** (near-uniform) |

Three of the five zonal parameters carry a concrete, citable discrepancy (kS, kGSH, bG); the other two (k450, kG) are confirmed and now have numeric magnitudes.

## Caveats

- DVP coefficient (protein regression) and Yakubovsky log2-ratio (mRNA) are on different scales; use one scale consistently when parameterizing, don't mix. Prefer DVP (protein) for reaction rates.
- CYP3A4 and UGT2B7 are human-specific; they appear only in the DVP protein set, not in the cross-species mRNA comparison.
- Values are *baseline healthy-liver* zonation. None of these datasets measure APAP-exposed or diseased zonation, so recalibration applies baseline gradients to the existing kinetic core (see `spatial_liver_zonation_datasets.md` open question #3).

---

## CORRECTION (2026-09-13, reviewer-mandated)

The DVP "coefficient" is β1 from a robust-scaled linear mixed model (median/IQR
scaling), **not** a portal/central abundance ratio — the reviewer is right. The
same table's 20-bin intensity profiles (which sum to ~100%) *do* give an
abundance ratio. Re-derived folds (central-end / portal-end bins, protein-level):

| Gene | robust fold C/P | Gene | robust fold C/P |
|---|---|---|---|
| CYP1A2 | 6.9× | UGT1A1 | 2.1× |
| CYP3A4 | 3.5× | UGT1A6 | 1.8× |
| CYP2E1 | **1.9×** | UGT2B7 | **1.4×** |
| SULT1A1 | 0.99× (uniform) | GSTA2 | **1.16×** |
| GCLC/GCLM/GSS | ~1.0× (uniform) | | |

Protein folds are systematically lower than the mRNA (Yakubovsky) log2 ratios
(e.g. GSTA2 1.16× protein vs 3.4× mRNA) — a genuine mRNA–protein discordance.
Protein abundance is itself a *proxy* for catalytic capacity (assumption).
Model now uses these protein folds with **arithmetic-mean**-preserving ramps
(total enzyme conserved), per the reviewer's point 5.

**Consequence for the results:** the dose-dependent sign flip did NOT survive —
recalibration is now *protective at both doses* (see `paper_outline.md` update).
