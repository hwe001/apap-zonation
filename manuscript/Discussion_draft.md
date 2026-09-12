## Discussion

**The stability boundary is the finding.** Replacing the assumed zonation
gradients with measured human profiles produces a sharply dose-dependent answer:
at the therapeutic dose (4 g) the predicted pericentral adduct burden falls by
~90% (90% CI −97 to −81%), while at overdose (16 g) the same replacement is null
(median −4%, CI −20 to +11%). The assumed gradients therefore overestimate
pericentral injury specifically at low dose, and are irrelevant at high dose.
This is not a recalibration result that happens to be fragile; it is a stability
result, obtained precisely because we propagated uncertainty rather than fixing a
single set of folds.

**Why the boundary falls where it does.** The two corrections that matter are the
glutathione-synthesis gradient and the CYP gradient. The assumed gradients place
glutathione synthesis strongly periportal (9×), so pericentral cells begin with
far less glutathione; measured GCLC/GCLM/GSS are essentially unzonated, so
pericentral glutathione is higher than assumed. The assumed CYP gradient (5×) is
also steeper than the measured CYP2E1 fold (1.9×), although CYP1A2 (6.9×) is
steeper, which is why the CYP correction contributes uncertainty rather than a
clean sign. At the therapeutic dose these two corrections reduce pericentral
NAPQI formation and increase its detoxification, cutting adducts ~90%. At
overdose, glutathione is depleted to near-zero in every cell regardless of its
zonation, so the synthesis-gradient correction becomes moot and the prediction is
governed by total NAPQI formation, which is insensitive to the gradient shape.
The one correction that pushes the other way — glucuronidation, measured ~1.4×
rather than the assumed 9× — is real but small relative to the glutathione effect
at low dose, and is why the overdose result is a null rather than a robust
reduction.

**Comparison with the state of the art.** Prior zonated acetaminophen models —
our own Means & Ho (2019) and Franiatte et al. (2019), and the 5,114-cell lobule
of Ghosh et al. (2026) — all impose the same qualitative, Gebhardt-derived
gradients and report point predictions from them; none has tested which of their
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
The proteomic folds are cohort-level; donor-level uncertainty was approximated by
the protein–mRNA span rather than measured directly, so the confidence intervals
in Fig. 4 are a lower bound on the true uncertainty. The model is a 16-cell
sinusoid without transport, so it quantifies the consequences of the gradients
but cannot resolve transport-coupled questions; we deliberately frame it as a
zonation-sensitivity study rather than a full reconstruction. The adduct
compartment is a proxy for injury, not necrosis itself, and its elimination rate
is a calibrated input, not a validated output.

**Conclusion.** For a model of acetaminophen hepatotoxicity, the answer to "how
much does assumed zonation matter?" is: a great deal at the therapeutic dose, and
almost not at all at overdose. The stability test that yields this answer is the
transferable contribution — a way to ask, of any zonated toxicity model, which
predictions are robust to the one assumption that has been least examined.
