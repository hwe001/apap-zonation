## Discussion

**The stability boundary is the finding.** Replacing the assumed zonation
gradients with measured human profiles produces a sharply dose-dependent answer:
at the simulated 4 g single dose the predicted pericentral adduct burden
falls by ~93% (90% simulation interval −96 to −89%), at the simulated 8 g single
dose by ~37% (interval −43 to −30%), while at the simulated 16 g single dose the
net response is small and uncertain (median −7%, interval −13 to +0.3%). The
assumed gradients therefore overestimate the pericentral adduct burden at low
dose, but at high dose their replacement produces large, offsetting effects
rather than a robust net change. This is not a recalibration result that happens
to be fragile; it is a stability result, obtained precisely because we propagated
uncertainty rather than fixing a single set of folds.

**Why the boundary falls where it does.** The two corrections that matter are the
glutathione-synthesis gradient and the CYP gradient. The assumed gradients place
glutathione synthesis strongly periportal (9×), so pericentral cells begin with
far less glutathione; measured GCLC/GCLM/GSS are essentially unzonated, so
pericentral glutathione is higher than assumed. The assumed CYP gradient (5×) is
also steeper than the isoform-weighted measured fold (3.4×); because CYP2E1
dominates human APAP oxidation, the abundance weighting narrows the CYP
uncertainty to a modest, consistently protective correction. At the simulated 4 g single dose these two corrections reduce pericentral
NAPQI formation and increase its detoxification, cutting adducts ~90%. At
overdose, glutathione is depleted in every cell, which weakens the
glutathione-synthesis correction (−27% at 16 g versus −83% at 4 g), so it no
longer dominates; instead it is nearly cancelled by the glucuronidation
correction, which is substantial (+50%) and pushes the other way because a
shallower glucuronidation gradient leaves more APAP available for oxidation. The
result is not that zonation ceases to matter at overdose — each gradient still
moves the prediction by tens of percent — but that the individual effects oppose
each other, so their sum is small and uncertain.

**Comparison with the state of the art.** Prior zonated acetaminophen models —
our own Means & Ho (2019) and Franiatte et al. (2019), which impose
Gebhardt-derived linear gradients, and the 5,114-cell lobule of Ghosh et al.
(2026), which builds zonation from hypothesis-based incremental percentages —
all construct their zonal enzyme profiles from the literature rather than
measurement, and report point predictions from them; none has tested which of
their
conclusions survive replacing the gradient with measurement, and the largest
states explicitly that no genomic or proteomic data were integrated. Whole-liver
quantitative-systems-toxicology platforms such as DILIsym omit zonation entirely,
which sidesteps rather than answers whether the spatial detail matters. Our
contribution is therefore twofold. The *empirical* novelty is that measured human
zonation is systematically shallower than the imposed gradients and materially
changes low-dose predictions — a fact that any future zonated model must
confront. The *methodological* novelty is the stability test itself: a way to
separate the robust from the assumption-dependent conclusions of a zonated
toxicity model, which we argue should accompany any model whose predictions
depend on zonal enzyme assumptions.

**Limitations.** The folds are protein abundances, and abundance is a proxy for
catalytic capacity, not capacity itself — no lobule-resolved activity data exist.
The proteomic folds are cohort-level; donor-level variation was not measured and
is not represented in the simulation intervals, so the intervals in Fig. 4
describe scenario uncertainty over fold ranges, not donor-level sampling
uncertainty. The model is a 16-cell
sinusoid without transport, so it quantifies the consequences of the gradients
but cannot resolve transport-coupled questions; we deliberately frame it as a
zonation-sensitivity study rather than a full reconstruction, and the stability
boundary established here holds for this reduced system — whether it persists in
transport-coupled lobule models is untested and is precisely what the stability
test is designed to assess. The adduct
compartment is a proxy for injury, not necrosis itself, and its elimination rate
is a calibrated input, not a validated output.

**Conclusion.** For a model of acetaminophen hepatotoxicity, the answer to "how
much does assumed zonation matter?" is dose-dependent and more subtle than a
single statement. At the simulated 4 g single dose the assumed gradients
substantially overestimate the pericentral adduct burden, and replacing them with
measured profiles reduces it by ~90%. At overdose (16 g) the individual gradient
corrections remain large — the glutathione-synthesis and CYP corrections reduce
adducts, the glucuronidation correction increases them — but they oppose each
other and nearly cancel, leaving a small, uncertain net response. The stability
test that yields this answer is the transferable contribution — a way to ask, of
any zonated toxicity model, which predictions are robust to the one assumption
that has been least examined.

## Code and data availability

The model code, the fold-derivation inputs, the figures, and an interactive
browser viewer are publicly available at
[github.com/hwe001/apap-zonation](https://github.com/hwe001/apap-zonation).
The submission version analysed in this paper is permanently available as
release v1.1.0
([github.com/hwe001/apap-zonation/releases/tag/v1.1.0](https://github.com/hwe001/apap-zonation/releases/tag/v1.1.0));
a DOI will be minted for this release via Zenodo on acceptance. The viewer
(`viewer.html`) displays precomputed simulations of the five model
species across the 16 hepatocytes, both doses, and both gradient schemes; it is
intended for inspecting the results and does not constitute validation evidence.
