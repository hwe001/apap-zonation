# JTB submission checklist

Status against the Journal of Theoretical Biology requirements and the reviewer's
second-round comments. "Done" means the item is in the current draft; "needs" flags
what is still required before submission.

## Reviewer round-2 items

| # | Item | Status |
|---|---|---|
| 1 | Relabel "90% CI" as a model-scenario/simulation interval | ✅ "90% simulation interval" throughout; fold ranges now stated explicitly in Methods §2.4 |
| 2 | Overdose finding reframed (opposing effects, not "irrelevant") | ✅ carried through Highlights, abstract, §3.3, Discussion, Conclusion |
| 3 | Fold derivation reproducible (bins, averaging, inputs) | ✅ `data/fold_derivation.md` + Methods §2.3 (bins 0–4 vs 15–19, mean-then-ratio) |
| 4 | Tone down injury claims; 4 g as "single dose" | ✅ "adduct burden" replaces "injury" in claims; "therapeutic single dose (4 g)" |
| 5 | Absolute predictions alongside % changes | ✅ §3.2 and §3.3 report absolute nmol and %-of-dose |

## JTB front matter

| Item | Status |
|---|---|
| Title | ✅ (stability question) |
| Highlights (biological application + theoretical advance) | ✅ in `manuscript/front_matter.md`; verify each ≤ 85 characters |
| Abstract (single paragraph) | ✅ |
| Keywords | ✅ |
| Code and data availability statement | ✅ in README + to be added to the manuscript |

## Manuscript body

| Item | Status |
|---|---|
| Introduction (novelty + state of the art) | ✅ |
| Methods (model, fold derivation, uncertainty) | ✅ |
| Results (stability finding, absolute + relative) | ✅ |
| Discussion (novelty, comparison, limitations) | ✅ |
| Conclusion | ✅ |
| References (complete, no markers) | ✅ (Means & Ho completed as 2019) |
| Figures 1–4 (self-contained) | ✅ in `figures/` |
| Tables 1–3 | ✅ |

## Cleanup

| Item | Status |
|---|---|
| Remove working-draft notes | ✅ |
| Remove unresolved citation markers | ✅ |
| Remove internal filesystem paths | ✅ (internal notes gitignored) |

## Before submission — needs attention

- [ ] **Archive a versioned DOI** (e.g. Zenodo) of the repository and cite that
      version in the paper; keep the GitHub link pointing at the current version.
- [ ] **Render the bundled document** and check figure legibility / page layout —
      the reviewer could not render it in their environment.
- [ ] **Confirm Highlights ≤ 85 characters each** (JTB/Elsevier limit).
- [ ] **Set the LICENSE copyright holder** if "The Authors" is not final (done:
      Harvey Ho).
- [ ] **Cover letter** stating the theoretical advance and biological insight.
- [ ] **Abstract word count** within the journal limit.
- [ ] **Author list / affiliations / funding / conflict-of-interest** fields.
