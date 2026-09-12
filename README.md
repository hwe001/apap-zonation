# Acetaminophen zonation — a stability analysis

Companion repository for

> **How much does assumed enzyme zonation matter for predicted acetaminophen
> hepatotoxicity? A stability analysis against measured human profiles.**

This is a research prototype accompanying a manuscript under review. It is **not**
a validated clinical or drug-safety tool.

## What this repository shows

Acetaminophen (paracetamol) injury is spatially restricted to the pericentral
zone of the liver because the enzymes of its metabolism are zonated along the
hepatic sinusoid. Prior zonated models impose *literature-derived* gradients for
those enzymes; this work replaces them with *measured* human single-cell spatial
proteomics and asks a stability question: **which predictions survive the
replacement, and which do not?**

The headline result is dose-dependent: at the therapeutic single dose (4 g) the
assumed gradients overestimate the pericentral adduct burden (~90% over-correction),
while at overdose (16 g) the individual gradient corrections remain large but
**oppose each other** — the shallower glutathione-synthesis and CYP gradients
reduce adducts, the shallower glucuronidation gradient increases them — so the
net response is small and uncertain.

## Try it in your browser — no installation

Open **[`viewer.html`](viewer.html)** in any modern browser. It renders the
per-hepatocyte time courses as interactive heatmaps with:

- a **time slider / play** to animate the concentration wave along the sinusoid,
- a **species** selector (APAP, PAPS, NAPQI, glutathione, protein adducts),
- an **assumed vs measured** gradient toggle,
- a **dose** selector (4 g and 16 g).

The data are precomputed and embedded, so the file is fully self-contained. It is
a companion for *inspecting* the results, not evidence that validates the model —
the manuscript's figures and numerical results are self-contained.

## Reproduce the results

```bash
pip install -r requirements.txt
python code/zonated_apap_model.py   # model + one-at-a-time decomposition
python code/uncertainty.py          # Monte Carlo uncertainty (saves results/*.npz)
python code/make_figures.py         # regenerate figures/ (needs uncertainty.py first)
python code/make_viewer.py          # regenerate viewer.html
```

## Repository structure

```
viewer.html                    interactive browser viewer (self-contained)
code/
  intracellular_apap_model.py  single-hepatocyte APAP kinetics (Reddyhoff 2015)
  zonated_apap_model.py        16-hepatocyte zonated model + gradient schemes
  uncertainty.py               Monte Carlo fold-uncertainty propagation
  make_figures.py              publication figures
  make_viewer.py               viewer generation
figures/                       the four manuscript figures (PNG)
manuscript/                    manuscript draft (markdown)
data/                          measured zonation coefficients
```

## Data provenance

The measured zonation folds are re-derived from the binned intensity profiles of
the single-cell deep-visual proteomics dataset (Weiss et al., *Nat Metab* 2026),
not from its published regression coefficients (which are not fold changes).
See `data/measured_zonation_coefficients.md` for the conversion and the full
table. Protein abundance is a proxy for catalytic capacity; no lobule-resolved
activity data exist.

## Code and data availability

The model code, fold-derivation inputs, figures, and the interactive viewer are in
this repository. `viewer.html` displays **precomputed** simulations (not live
model runs) of the five species across the 16 hepatocytes, both doses, and both
gradient schemes. The measured zonation folds are re-derived from the
single-cell deep-visual proteomics bin profiles (Weiss et al., *Nat Metab* 2026)
as documented in `manuscript/Supplementary_Methods.md` (per-protein bin values,
calculated ratios, and uncertainty inputs) and summarised in
`data/measured_zonation_coefficients.md`. A versioned archive (DOI) of the
submission version will be cited in the paper; the GitHub link points to the
current version.

## License

MIT. See `LICENSE`.

## Citation

If you use this work, please cite the accompanying manuscript (citation to be
added on publication).
