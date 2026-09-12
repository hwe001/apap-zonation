## Introduction

Acetaminophen (paracetamol) is safe at therapeutic doses but a leading cause of
drug-induced acute liver failure in overdose, through a minor oxidative pathway
that forms the reactive metabolite *N*-acetyl-*p*-benzoquinone imine (NAPQI),
which depletes glutathione and covalently binds hepatocyte proteins (Reddyhoff
et al. 2015). The injury is spatially restricted to the pericentral (zone 3)
region because the enzymes involved are zonated along the periportal–pericentral
axis: the cytochrome P450 isoforms that generate NAPQI are expressed
predominantly pericentrally, while sulfation, glucuronidation and glutathione
conjugation have their own zonal patterns (Gebhardt 1992; Halpern et al. 2017).

Because this zonation determines where and when injury occurs, it has been built
into a sequence of mathematical models of increasing complexity. Means & Ho
(2019) embedded a five-species intracellular kinetic model (Reddyhoff et al.
2015) in a finite-element sinusoid; Franiatte et al. (2019) added
advective-diffusive blood-flow transport; and a separate group has since built a
full three-dimensional lobule of 5,114 hepatocytes coupled to a whole-body PBPK
model (Ghosh et al. 2026). All three, however, share one assumption: the zonal
gradients for the enzyme rates are *imposed* as qualitative, literature-derived
profiles, ultimately traceable to a generic, largely rodent-derived review
(Gebhardt 1992), rather than measured. The largest current model states this
explicitly — "no genomic or proteomic data were directly integrated" (Ghosh et
al. 2026) — and the earliest was criticised for exactly this at review.

Two 2026 human spatial-omics datasets now make the gradients measurable:
single-cell deep-visual proteomics (Weiss et al. 2026) and a multi-modal spatial
transcriptomics atlas (Yakubovsky et al. 2026). This raises a question that is at
once empirical and, we argue, theoretical. The straightforward contribution
would be to "recalibrate" the gradients and report the resulting prediction; but
a single recalibrated prediction inherits the fragility of whichever fold values
and isoform weightings are chosen, and does not tell a biologist which prior
conclusions still hold. We therefore ask a **stability question**: which
conclusions about pericentral acetaminophen injury survive replacing the assumed
gradients with measured human profiles, and which do not?

We answer it in three steps. First, we derive per-enzyme abundance folds from the
proteomic bin profiles — not from the published regression coefficients, which
are not fold changes — and document the conversion and its uncertainty (Section
2.2). Second, we impose these folds on the 16-hepatocyte model with total enzyme
conserved, and propagate isoform and protein–mRNA uncertainty by Monte Carlo
(Sections 2.3–2.5). Third, we report which predictions survive (Section 3). The
result is a dose-dependent stability boundary: at the therapeutic single dose the
assumed gradients robustly overestimate the pericentral adduct burden, whereas at
overdose the individual gradient corrections remain large but oppose each other,
so the net response is small and uncertain. We close by discussing what this
stability test contributes as a method for zonated toxicity models generally
(Section 4).
