## Methods

### 2.1 Model

The per-hepatocyte kinetics follow the five-species acetaminophen (APAP)
metabolism model of Reddyhoff et al. (2015), which tracks paracetamol *P*, the
sulfation co-substrate PAPS *S*, the reactive metabolite NAPQI *N*, glutathione
*G*, and covalent protein adducts *C*:

```
dP/dt = −kS·S·P − kG·P − k450·P + kN·N
dS/dt = −kS·S·P + bS − dS·S
dN/dt = k450·P − kN·N − kGSH·N·G − kPSH·N
dG/dt = −kGSH·N·G + bG − dG·G
dC/dt = kPSH·N − k_clear·C
```

The first five terms are those of Reddyhoff et al.; we add a first-order
elimination term to the adduct compartment so it has the correct turnover. Its
rate `k_clear = 0.42 d⁻¹` is the measured human serum adduct-elimination rate
(James et al. 2009); we stress that this is a *calibrated input*, not an
independent test, and it enters only the adduct turnover, not the zonal analysis
that is the subject of this paper.

Parameter values were taken from the original MATLAB source rather than the
published table (whose parameter names appear transposed relative to their
values), and reproduce the paper's stated dosing conversion exactly (Table 1).
The model is human-scaled: dosing (4 g therapeutic, 16 g overdose), hepatocyte
number (1.6055 × 10¹¹), and the glucuronidation/oxidation rates are human-derived
(Reddyhoff et al. 2015).

**Table 1 — Kinetic parameters (Reddyhoff et al. 2015).**

| Parameter | Value | Meaning |
|---|---|---|
| kS | 2.26 × 10¹⁴ · 10⁻¹² cell·mol⁻¹·d⁻¹ | PAPS-dependent sulfation |
| kG | 2.99 d⁻¹ | glucuronidation |
| k450 | 0.315 d⁻¹ | CYP450 oxidation |
| kN | 0.0315 d⁻¹ | NAPQI → APAP reverse |
| bS | 2.65 × 10⁻¹⁴ / 10⁻¹² mol·cell⁻¹·d⁻¹ | PAPS synthesis |
| dS | 2.0 d⁻¹ | PAPS decay |
| kGSH | 1.6 × 10¹⁸ · 10⁻¹² cell·mol⁻¹·d⁻¹ | GSH conjugation of NAPQI |
| kPSH | 110 d⁻¹ | NAPQI–protein binding |
| bG | 1.374 × 10⁻¹⁴ / 10⁻¹² mol·cell⁻¹·d⁻¹ | GSH synthesis |
| dG | 2.0 d⁻¹ | GSH decay |
| k_clear | 0.42 d⁻¹ | adduct elimination (James 2009) |

### 2.2 Zonated extension

The intracellular model is run independently across *N* = 16 hepatocytes arrayed
along a periportal→pericentral sinusoid (hepatocyte 1 periportal, 16
pericentral). Each hepatocyte receives the same initial APAP dose
`P0 = dose · f / (MW · Nhep)`, with hepatic penetration `f = 0.8`, molecular
weight `MW = 151 g·mol⁻¹`, and `Nhep = 1.6055 × 10¹¹` hepatocytes. Enzyme rates
are modulated by a linear ramp along the sinusoid. For an enzyme with
pericentral/periportal fold *F*, the per-hepatocyte multiplier is

```
m(x) = a + (b − a)·x,   x = (i − 1)/15,   a = 2/(1+F), b = 2F/(1+F),
```

so that the *arithmetic* mean multiplier across the sinusoid is exactly 1 — total
enzyme is conserved and only its location varies. (A geometric-mean-preserving
ramp would not conserve total enzyme; we use the arithmetic form throughout.)

[[FIG fig1_schematic.png | Model schematic. (Left) Sixteen hepatocytes arrayed along a periportal-to-pericentral sinusoid, each running the same intracellular kinetics with enzyme rates modulated by a zonal gradient. (Right) The five-species APAP network: APAP (P) is cleared by sulfation (SULT) and glucuronidation (UGT), or oxidized by CYP450 to NAPQI (N), which is detoxified by GSH conjugation (GST) or binds protein to form the toxic adducts (C).]]

### 2.3 Measured zonation and fold derivation

We used the 2026 single-cell deep-visual proteomics (DVP) dataset of human liver
zonation (Weiss et al. 2026). Two features of that dataset matter for the
conversion, and both are worth stating plainly. First, the published "zonation
coefficient" is the fixed-effect slope of a per-protein linear mixed model fitted
to robust-scaled (median-centred, IQR-divided) intensities; it is *not* a
portal-to-central abundance ratio, and cannot be read as an enzyme-rate fold
change. Second, the same supplementary table reports each protein's mean
intensity in 20 equidistant bins along the porto-central axis, normalised so the
bins sum to 100%. We therefore re-derived abundance folds directly from these bin
profiles: for each protein we averaged the five central-most bins (bins 0–4) and
the five portal-most bins (bins 15–19), then took the central/portal ratio of
those two averages (Table 2; full per-protein inputs in the Supplementary
Methods). The DVP pipeline reports only proteins with ≥70% data completeness, and
the bin profiles are cohort-level means across the 14 healthy donors, so
donor-level variation is not propagated at this stage — it is approximated in the
uncertainty analysis (Section 2.4). Protein folds were cross-checked against mRNA log2(portal/central)
ratios from the Yakubovsky et al. (2026) atlas; the two agree in direction but
differ in magnitude (a genuine mRNA–protein discordance), which we carry into the
uncertainty analysis. We note explicitly that protein abundance is a *proxy* for
catalytic capacity, not capacity itself; no lobule-resolved activity data exist.

**Table 2 — Assumed versus measured zonation (pericentral/periportal fold).**

| Parameter | Enzyme | Assumed (prior models) | Measured (protein bin-ratio) |
|---|---|---|---|
| k450 | CYP450 | 5× pericentral | 3.4× (isoform-weighted: CYP2E1/1A2/3A4 = 0.60/0.25/0.15; raw isoform folds 1.9–6.9×) |
| kG | UGT | 9× pericentral | 1.4× (UGT2B7; range 1.4–2.1×) |
| kS | SULT | 9× periportal | unzonated (SULT1A1 0.99×) |
| kGSH | GST | 9× periportal | 1.2× pericentral (GSTA2; reversed direction) |
| bG | GSH synthesis | 9× periportal | unzonated (GCLC/GCLM/GSS ~1.0×) |

### 2.4 Gradient replacement and stability analysis

The five enzyme parameters (k450, kG, kS, kGSH, bG) were assigned the measured
folds; sulfation (kS) and GSH synthesis (bG) are unzonated, and GSH conjugation
(kGSH) is weakly pericentral rather than periportal. Because the CYP isoforms
contribute very unequally to human APAP oxidation — CYP2E1 dominates, CYP1A2 is
secondary, CYP3A4 minor (Laine et al. 2009) — the k450 fold is the
isoform-weighted combination 0.60·(CYP2E1 fold) + 0.25·(CYP1A2) + 0.15·(CYP3A4)
= 3.4×, rather than the raw isoform span. To test robustness, we propagated the
fold uncertainty by 300-draw Monte Carlo at three simulated doses (4, 8 and
16 g). The k450 fold was sampled by drawing the isoform weights over plausible
ranges (CYP2E1 0.50–0.70, CYP1A2 0.15–0.30, CYP3A4 the renormalised remainder)
and recomputing the weighted fold; kG was sampled log-uniformly over [1.4, 2.1]
(UGT2B7-to-UGT1A1); kGSH over [1.16, 3.4] (the protein-to-mRNA span for GSTA2);
and kS, bG over [0.9, 1.1] (±10% around uniform — a stated assumption for the
width, not an estimated dispersion). We report the
distribution of the resulting change in pericentral adduct burden relative to the
assumed gradients. The assumed
(literature) gradients are held fixed as the null hypothesis being tested; only
the measured folds are varied. Because the 5×/9× assumed folds are themselves
stylised (the prior models specify qualitative "strong" gradients), we also
re-ran the 4 g comparison against 27 variants of the assumed gradients (CYP
3–7×; UGT and the periportal gradients 6–12×).

### 2.5 Numerical analysis

The coupled ODE system — 16 hepatocytes × 5 intracellular species = 80 state
variables — was integrated with a stiff Radau solver
(`solve_ivp`, rtol 10⁻⁵, atol 10⁻⁸) over 3 days (the adduct peak occurs within the
first day; results are unchanged at tighter tolerances). Reported quantities are
the peak adduct concentration per hepatocyte (C_pk) and its pericentral value.
The one-at-a-time attribution moves each parameter from its assumed to its
measured gradient while holding the others fixed; the Monte Carlo varies all
folds jointly.

---

## Results

### 3.1 Measured zonation is mostly shallower than assumed

The measured human zonation disagrees with the literature-derived gradients in
direction and magnitude (Table 2, Fig. 2). Three enzymes are unzonated or nearly
so — sulfation (SULT1A1, fold 0.99×) and glutathione synthesis (GCLC/GCLM/GSS,
folds ~1.0×) — where the prior models imposed strong periportal gradients. The
glutathione-conjugation gradient is reversed in direction (GSTA2 is weakly
pericentral, 1.2×, rather than strongly periportal). The two gradients that keep
their assumed direction — CYP450 and glucuronidation — are both shallower than
assumed: weighted by the isoforms' approximate contributions to human APAP
oxidation (CYP2E1 ≈ 0.60, CYP1A2 ≈ 0.25, CYP3A4 ≈ 0.15; Laine et al. 2009), the
CYP fold is 3.4× versus the assumed 5×, and glucuronidation 1.4× versus 9×.
CYP1A2's raw protein fold (6.9×) does exceed the assumed 5×, but its minor
contribution to APAP oxidation leaves the weighted CYP fold below the
assumption; isoform uncertainty therefore enters through the weighting (Section
2.4) rather than through the raw isoform span.

[[FIG fig2_gradients.png | Assumed versus measured zonation gradients for the five enzyme parameters, as a fold-change relative to a uniform baseline (1×). Orange dashed: the literature-derived gradient carried by prior models; blue solid: the gradient re-derived from the 2026 proteomic bin profiles. Sulfation and GSH synthesis are unzonated; GSH conjugation is reversed (weakly pericentral); CYP450 and glucuronidation retain direction but are far shallower than assumed.]]

### 3.2 Gradient replacement reduces predicted pericentral adducts at the simulated 4 g single dose

Replacing the assumed gradients with the measured ones reduces the predicted
pericentral adduct burden at the simulated 4 g single dose by 90% (Table 3):
the peak pericentral adduct falls from 1.1 × 10⁻³ to 1.1 × 10⁻⁴ nmol per
hepatocyte (0.83% to 0.085% of the dose). The
reduction is dominated by flattening the glutathione-synthesis
gradient (−83%), which raises pericentral glutathione; the CYP
correction (−8.5%) and the sulfation correction (−9%) are secondary, and the
glucuronidation
correction (+52%) acts in the opposite direction — a shallower glucuronidation
gradient leaves more APAP available for oxidation — but is outweighed by the
glutathione effect.

[[FIG fig3_damage.png | Predicted peak adduct fraction of dose along the sinusoid at 4 g and 16 g, for the assumed (orange dashed) and measured (blue solid) gradients. At the simulated 4 g single dose the measured gradients predict a substantially lower, less sharply zoned pericentral adduct burden; at the simulated 16 g single dose the two are similar.]]

### 3.3 At overdose, large opposing corrections nearly cancel

The correction is strongly dose-dependent. At the simulated 4 g single dose
the measured gradients reduce the pericentral adduct burden by a median 93%
(90% simulation interval −96 to −89%), and this is robust: all 300 Monte Carlo
draws give a reduction. An intermediate simulated dose (8 g) shows a smooth
transition (median −37%, interval −43 to −30%; 100% protective). At the
simulated 16 g single dose the net response is small and uncertain (median −7%,
interval −13 to +0.3%), with 94% of draws protective and 6% aggravating
(Fig. 4); the point-estimate peak pericentral adduct moves from 5.5 × 10⁻³ to
5.6 × 10⁻³ nmol per hepatocyte (1.04% to 1.06% of the dose). This is not because
zonation ceases to matter at overdose: each gradient still moves the prediction
substantially (Table 3). Rather, the protective corrections — the
glutathione-synthesis gradient (−27%) and, to a lesser extent, the CYP gradient
(−7%) — are nearly
cancelled by the aggravating glucuronidation correction (+50%), because a
shallower glucuronidation gradient leaves more APAP available for oxidation. The
net response is small and uncertain precisely because large effects oppose each
other. The conclusion is also insensitive to the exact stylised values of the
assumed gradients: re-running the 4 g comparison against 27 variants of the
assumed folds (CYP 3–7×; UGT and the periportal gradients 6–12×) gives a median
reduction of −89% (range −92 to −86%), protective in all variants.

[[FIG fig4_uncertainty.png | Stability of the gradient replacement under uncertainty (300 Monte Carlo draws per dose over the isoform-weighted CYP, UGT and protein–mRNA fold ranges). Left to right: simulated 4, 8 and 16 g single doses. The reduction is robust at 4 g (median −93%, 90% simulation interval −96 to −89%) and 8 g (−37%, −43 to −30%; 100% of draws protective in both), and fades to a small, uncertain net response at 16 g (−7%, −13 to +0.3%; 94% protective / 6% aggravating).]]

**Table 3 — One-at-a-time effect on pericentral peak adducts (simulated single doses).**

| Change (assumed → measured) | Δ at 4 g | Δ at 8 g | Δ at 16 g |
|---|---|---|---|
| kG: glucuronidation 9× → 1.4× | +52% | +50% | +50% |
| bG: GSH synthesis periportal → unzonated | −83% | −50% | −27% |
| k450: CYP 5× → 3.4× (isoform-weighted) | −8.5% | −7.6% | −7.2% |
| kS: sulfation periportal → unzonated | −9% | −4% | −1% |
| kGSH: GST flip → weakly pericentral | −4% | −2% | −1% |
