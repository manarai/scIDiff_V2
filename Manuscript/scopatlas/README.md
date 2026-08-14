# scOpAtlas manuscript — Bioinformatics Applications Note draft

Reframed methods/software note after the reviewer critique of the
prior four-regime biology-paper draft. This version:

- Corrects the tercile-classifier scientific error (absolute
  hyperbolicity primary + data-driven secondary axes; tercile retained
  only as `thresholds='rank'` opt-in).
- Includes a **replicated empirical observation** across TWO
  hematopoiesis datasets (Paul15 + Setty 2019 marrow HSPC): ≥97% of
  cells are non-hyperbolic under standard scJDO training. This is the
  paper's honest headline; two candidate explanations (biology vs
  velocity-prior bias) are flagged as open.
- Adds a low-drift-artifact diagnostic (Spearman(P, ‖f‖) = −0.48
  Paul15, −0.31 Setty) that reframes the *plastic* regime as a
  low-drift-region indicator rather than a biological plasticity
  claim.
- Reframes ensemble concordance honestly (18.7% discordant on Paul15
  with 5-member ensemble varying seed × velocity-prior).

## Contents

- **`manuscript.md`** — App-Note main text (~2,700 words).
  Abstract → intro → methods → results (3 subsections) → discussion →
  availability → 1 main figure caption → 6 supplementary figures
  captions.
- **`references.md`** — bibliography with placeholders where DOIs need
  verification.
- **`figures/`**:
  - `fig1_appnote.png` — **main figure**: 3-panel headline.
    (A) λ_max⁺ histograms Paul15 vs Setty (the ≥97% observation);
    (B) plasticity vs drift magnitude on both (the low-drift artifact);
    (C) Paul15 ensemble concordance (19% discordant at 0.5).
  - `figS1_synthetic_validation.png` — 3-panel synthetic validation
    (regime map + concordance + negative control).
  - `figS2_paul15_regime_umap.png` — Paul15 UMAP colored by regime
    (rank mode) + pseudotime reference.
  - `figS3_paul15_markers.png` — marker overlay + regime overlay on
    Paul15 UMAP. **Caveat**: overlaps with the low-drift artifact.
  - `figS4_paul15_regime_by_stage.png` — regime × stage per lineage
    (Paul15, rank mode). Within-model consistency check.
  - `figS5_paul15_threshold_sensitivity.png` — regime-fraction
    sensitivity to (τ_upper, deep_cut) grid.
  - `figS6_paul15_concordant_vs_all.png` — non-redundancy comparison.
  - `raw/` — original individual panel PNGs, older 6-panel Fig 1
    composite, plus the standalone synthetic sub-panels used to build
    Fig S1.

## Before you submit

1. **Author list is `[Author A]` / `[Author B]` / `[Corresponding
   Author]`**. Do NOT submit without confirming authorship with the
   scJDO team.
2. **Verify all DOIs in `references.md`.** Some are placeholder-quality.
3. **The two-dataset observation is the paper's headline finding.**
   Before final submission, run an ablation across `vel_scale ∈ {0,
   0.5, 2, 5}` on both datasets and report the per-cell fraction of
   hyperbolic cells as a function of `vel_scale`. This distinguishes
   the two candidate explanations (real biology vs training bias) and
   converts the "open question" into a resolved answer. Estimated
   effort: ~2 hours of compute + a supp figure.
4. **Consider adding a marrow-terminal dataset** (fully-differentiated
   tissue where hyperbolicity SHOULD exist a priori) to show the
   classifier does discriminate when biology permits. Would materially
   strengthen the App-Note.
5. **Regenerate figures at 300 DPI** for camera-ready submission.
6. **Author response note**: reviewer critique of the prior draft
   (auto tercile, plastic/deeply-stable overlap, missing depth-matched
   test, missing SpliceJAC/scMomentum, ensemble conflates gauge with
   spec) is addressed as follows: (1) classifier default changed to
   absolute hyperbolic; (2) plastic/deeply-stable now demonstrated as
   overlapping low-drift regions (Fig 1B) rather than distinct
   biological states; (3) depth-matched test is available in library
   but explicitly NOT claimed to be validated on Paul15 (§3.4); (4)
   SpliceJAC/scMomentum cited in §1; (5) gauge-vs-spec conflation
   acknowledged in §4 with rigorous per-cell gauge ensemble flagged
   as future work.

## Reproducibility

- `examples/scopatlas/tutorial.ipynb` — synthetic validation.
- `examples/scopatlas/paul15_test.ipynb` — Paul15 real-data workflow.
- `tests/test_scopatlas.py` — 48 tests including 2 CI acceptance gates
  (`test_ACCEPTANCE_negative_control`,
  `test_ACCEPTANCE_untrained_drift_fails_gate`).
- Setty marrow analysis was ad-hoc for this manuscript
  (`/scratchpad/marrow_investigate.py`) and should be promoted to
  `examples/scopatlas/setty_marrow_test.ipynb` before submission.

## Target venue

**Bioinformatics Applications Note** (primary). Alternatives: JOSS
(shorter, software-only, no App-Note-style biology claims), Genome
Biology Method (if a companion biology paper is added to the same
submission).
