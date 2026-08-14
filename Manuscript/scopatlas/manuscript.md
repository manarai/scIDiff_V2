# scOpAtlas: an operator-regime toolbox for single-cell drift fields, with dataset-scale calibration and gauge-uncertainty quantification

**[Author A]¹, [Author B]¹, [Corresponding Author]¹\***

¹[Affiliation]  ·  \*correspondence: [email]

*Placeholder author list. Do not submit without confirming authorship.*

---

## Abstract

Single-cell trajectory-inference methods (RNA velocity, CellRank, Palantir, dynamo, scJDO) infer a drift field $f(x,t)$ whose local Jacobian $J = \partial f/\partial x$ characterizes state-space stability. We present **scOpAtlas**, an open-source Python module (`scjdo.atlas`) that classifies cells into operator regimes from the per-cell Jacobian spectrum. The library is built around three engineering decisions: **(i) an absolute-hyperbolicity primary threshold** ($\lambda_{\max}^+ > 0$) that keeps regime labels comparable across datasets — replacing the tercile-based thresholds prior implementations use, which pin the unstable fraction at 33% by construction and destroy cross-cohort comparability; **(ii) ensemble concordance** across drift-field members varying random seed and velocity-prior configuration, written to `adata.obs['regime_concordance']` as a per-cell uncertainty measure; **(iii) a drift-quality gate** with three criteria (transport error, velocity-prior agreement, ensemble concordance) that runs before classification and blocks silent failure on under-fit models. Applied to two canonical mouse hematopoiesis datasets (Paul15, Setty 2019 HSPC), we report a replicated empirical observation: standard scJDO training yields ≥97% of cells with `λ_max⁺ > 0` — no fully-hyperbolic cells. We provide the library, tests, functional API, and reproducible notebooks; users should inspect their trained drift's spectrum before assuming a four-regime partition is meaningful on their data.

**Availability:** `github.com/manarai/scJDO` · `pip install scjdo` · MIT · Python ≥ 3.9.

---

## 1. Introduction

Trajectory-inference tools infer a drift field $f_\theta$ from which the local Jacobian $J(x,t)$ can be computed, and its eigenvalue spectrum characterizes cell-state stability, plasticity, and bifurcation proximity. Prior implementations that assign cells to dynamical regimes (dynamo [1], SpliceJAC [2], scMomentum [3]) use hand-set thresholds on eigenvalue-derived metrics. Three issues limit their reusability:

**(i) Threshold miscalibration.** Fixed thresholds do not transfer across ambient dimensions — order statistics of independent samples inflate with $n$. But data-driven thresholds tempt an obvious wrong fix: setting $\tau_{\text{upper}} = \text{quantile}_{2/3}(\lambda_{\max}^+)$ pins the unstable fraction at exactly 33% *by construction*, so a fully-hyperbolic dataset still gets a third of cells labelled unstable. Regime labels then cease to be comparable across datasets or conditions.

**(ii) Gauge freedom.** Snapshot single-cell data determines $f$ only up to a divergence-free circulation [4, 5]. Practical trainers resolve this by minimizing kinetic energy — a regularization choice, not a measurement. Individual complex eigenvalue pairs and precise bifurcation locations are gauge-variant; the SIGN of $\text{Re}(\lambda_{\max})$ is topologically robust (loss-of-hyperbolicity is invariant under smooth deformations). A regime classifier that reports point estimates without acknowledging this confuses gauge ambiguity with biological structure.

**(iii) Depth-confounded cohort tests.** The claim that regimes carry information beyond expression labels is testable as a regime × condition contingency effect within cell type. But regime depends on the embedding, which depends on UMI depth; without depth-matched controls the effect is trivially confounded.

We present **scOpAtlas** to address these three points at the library level: an absolute-hyperbolicity primary classifier that keeps labels dataset-comparable; an ensemble-concordance uncertainty quantifier; and a depth-matched permutation test with quality gate. The empirical section reports a replicated observation about scJDO-trained drifts on mouse hematopoiesis that users of the library should be aware of before interpreting regime output.

## 2. Methods

Given a trained drift field (any `torch.nn.Module` returning $f_\theta(x, t)$), scOpAtlas computes the per-cell Jacobian via `torch.func.jacrev` under `vmap` (~80× faster than row-by-row autograd; automatic fallback to classical autograd when the model contains ops that break batched tracing, e.g. kNN velocity lookups). Per cell we compute $\lambda_{\max}^+ = \max_j \text{Re}(\lambda_j)$, a stability quantile $q_\alpha = \text{quantile}_\alpha(\text{Re}(\lambda_j))$ (default $\alpha = 0.05$, dimension-stable proxy for $\lambda_{\min}^-$), a plasticity index $P = |\{j: |\lambda_j| < \varepsilon\}|/d$ (complex modulus), and the stable-subspace dimension $S = |\{j: \text{Re}(\lambda_j) < 0\}|$.

**Hyperbolic classifier (default `thresholds='hyperbolic'`).** Unstable is defined absolutely: $\lambda_{\max}^+ > 0$. Cells with all eigenvalues non-positive are classified further by data-driven scales tied to the median $|\text{Re}(\lambda)|$ across the dataset:

- **deeply-stable**: $q_\alpha < -k \cdot \text{median}|\text{Re}(\lambda)|$ (default $k = 2$)
- **plastic**: $|\lambda_{\max}^+| < k' \cdot \text{median}|\text{Re}(\lambda)|$ AND $P > 0.3$ (default $k' = 0.3$)
- **stable**: everything else

The tercile-based scheme is available as `thresholds='rank'` (kept for within-dataset ranking use cases) but is not default — it destroys cross-dataset comparability by construction.

**Ensemble concordance.** `OperatorEnsemble(drift_models, ...)` accepts $N \ge 2$ pre-trained drift models and writes per-cell concordance (fraction of members agreeing with the modal regime) to `adata.obs['regime_concordance']`. Members must vary the velocity-prior configuration (`vel_scale`, `vel_time_mode`) — seed-only variation samples a narrow slice of admissible drifts.

**Statistical inference.** `atlas.validate_nonredundancy(celltype_key, condition_key)` runs χ² + Cramér's V per cell type, permutation p-value from shuffling condition labels, Benjamini-Hochberg FDR across cell types, bootstrap CIs on regime fractions, and — critically — depth-matched subsampling (cells quantile-binned by UMI depth; each condition subsampled to equal count per bin; permutation repeated across many random matched draws). When `regime_concordance` is present, the entire pipeline runs twice — on all cells and on cells with concordance ≥ threshold — and both results are returned.

**Drift-quality gate.** Three criteria (transport error via pseudotime-holdout sliced-Wasserstein, RNA-velocity cosine agreement when available, ensemble concordance) run before regime classification. Each is skipped gracefully when its inputs are absent; overall pass requires ≥1 criterion ran and all that ran passed. On failure, `atlas.build()` raises `DriftQualityError`; `force=True` overrides with a loud banner.

**Test coverage.** 48 unit + integration tests; two acceptance gates pinned in CI (`test_ACCEPTANCE_negative_control` for pipeline calibration on synthetic data, `test_ACCEPTANCE_untrained_drift_fails_gate` for the quality-gate enforcement). Full API documented in `docs/scopatlas/`.

## 3. Results

### 3.1 Synthetic method validation (Supp. Fig. S1)

On 300 synthetic cells sampled along a 4D pseudotime trajectory with a hand-designed nonlinear drift, the classifier produces the intended regime partition, the 5-member ensemble concordance yields a non-trivial distribution with 14% discordant at threshold 0.5, and permutation inference under random condition assignment produces small non-significant Cramér's V values on both all-cells and concordant-only subsets — the acceptance-critical negative control for pipeline calibration. Identical-drift and wildly-different-drift ensemble cases are pinned in CI (`test_identical_members_give_perfect_concordance`, `test_ACCEPTANCE_negative_control`). Full walkthrough: `examples/scopatlas/tutorial.ipynb`.

### 3.2 Application to mouse hematopoiesis: two datasets, one replicated observation (Fig. 1A)

We applied scOpAtlas to two canonical mouse hematopoiesis scRNA-seq datasets: Paul15 [6] (2,321 cells across 15 clusters after removing low-count basophil/eosinophil/lymphoid populations) and Setty 2019 marrow HSPC [7] (4,142 cells, 11 Leiden clusters, DPT root at max-Cd34+ cell). For each dataset we trained a 5-member drift ensemble via `sjd.tl.fit_drift` (hidden=160, depth=3, 2000 epochs; members vary random seed and velocity-prior configuration).

Under the **hyperbolic classifier**, the per-cell $\lambda_{\max}^+$ distributions are essentially entirely positive:

| dataset | n | mean $\lambda_{\max}^+$ | 5%ile | 95%ile | fraction $\lambda_{\max}^+ < 0$ |
|---|---|---|---|---|---|
| Paul15 | 2,321 | +0.122 | +0.05 | +0.21 | **0.0%** |
| Setty marrow | 4,142 | +0.028 | +0.002 | +0.069 | **2.9%** |

**Neither dataset produces fully-hyperbolic cells under standard scJDO training.** All cells (Paul15) or 97% of cells (Setty marrow) have at least one direction of local instability. Fig. 1A shows the $\lambda_{\max}^+$ distributions overlaid; both bulk above zero.

This observation replicates across two independently trained drift fields on two datasets, so it is not a training-instance artifact. Two candidate explanations, both testable and both open:

1. **Real biology.** Hematopoietic differentiation is a continuous non-equilibrium flow — no cell sits at a Lyapunov-stable fixed point; every cell is transitioning somewhere. Terminal cells (e.g. mature erythrocytes) are in vivo constantly turning over. Under this reading, ≥97% non-hyperbolic is a genuine dynamical feature of hematopoiesis.

2. **Velocity-prior training bias.** scJDO's default velocity prior (`vel_scale=2.0`) biases the drift to follow the pseudotime gradient. Nonzero drift everywhere implies nonzero eigenvalues along the flow direction. A `vel_scale=0` ensemble member is included in our diversity design; regime call per cell is still dominated by the four velocity-primed members.

Distinguishing (1) from (2) requires either (a) a training ablation across `vel_scale ∈ {0, 0.5, 2, 5}` reported per-cell, or (b) a dataset where hyperbolicity should exist a priori (a fully-differentiated tissue snapshot). We flag this as an open question for the scJDO+scOpAtlas user community.

### 3.3 Plasticity captures low-drift regions, not near-neutral spectra (Fig. 1B)

We tested whether the *plastic* regime label — cells with $|λ_{\max}^+| < \varepsilon$ — captures a biological "near-neutral eigenvalue" state, or is a fitting artifact of low-drift regions where the network output is close to zero. Correlation on both datasets:

| dataset | Spearman(plasticity, $\|f(x)\|$) |
|---|---|
| Paul15 | −0.48 |
| Setty marrow | −0.31 |

Both datasets show a moderate-to-strong negative correlation between the plasticity index and local drift magnitude. Cells labelled *plastic* sit in low-drift regions of latent space. This is consistent with the hypothesis that the trained network output is near-zero in under-sampled regions and that plasticity indexes this rather than a distinct dynamical state. Users of the library should not interpret the `plastic` label as a biological plasticity claim without a separate control (e.g. matched non-neutral markers).

### 3.4 Ensemble concordance and non-redundancy (Fig. 1C, Supp. Figs S2–S6)

Under the rank-mode classifier retained for within-dataset ranking, Paul15 gives four populated regimes with per-cell ensemble concordance mean 0.64. **19% of cells fall below the 0.5 concordance threshold and 93% are not fully unanimous** across the 5-member ensemble — the honest gauge-uncertainty headline (Fig. 1C). Within-dataset non-redundancy on lineage × pseudotime-stage returns three significant lineages ($q < 0.05$) on Paul15, with a clean label-permutation negative control (0/5 nominal false positives under shuffled stage labels). Full per-cell UMAPs, per-lineage bar charts, marker validation cross-tables, and threshold-sensitivity sweeps in supplementary figures.

We note two caveats to reading these numbers as evidence of orthogonality to expression: (a) pseudotime and regime both derive from the same trajectory, so the "condition variable" is not independent of the model that produces the regime; (b) the ensemble diversity as constructed conflates gauge freedom with model-specification variance (`vel_scale=0` is a different estimator, not a gauge-admissible reparametrization). The depth-matched inference pipeline is available in the library (`validate_nonredundancy`) but is untested on Paul15 because the dataset has no cohort variable. A companion paper using a genuine cohort dataset is future work.

## 4. Discussion

scOpAtlas is a library, not a biological finding. Its scientific contribution — beyond the software — is a replicated empirical observation: canonical scJDO training on canonical mouse hematopoiesis produces essentially no fully-hyperbolic cells. This is a diagnostic finding for users. Before interpreting a four-regime partition as biological structure, users should:

1. Inspect the $\lambda_{\max}^+$ distribution and report the fraction < 0. If ≥95% of cells have $\lambda_{\max}^+ > 0$, the hyperbolic classifier collapses to a single label and the four-regime partition rests entirely on data-driven secondary thresholds. This is not necessarily wrong, but should be disclosed.
2. Compute Spearman(plasticity, $\|f(x)\|$). If strongly negative, the plastic regime is capturing low-drift regions rather than distinct spectral character. Rename or interpret accordingly.
3. Compute ARI(regime, cluster). A near-zero value indicates the classifier is not simply recapitulating expression labels; a substantial value indicates the two are entangled and downstream inference should account for this.

For gauge-uncertainty quantification, use ensemble concordance with members varying at minimum the velocity-prior configuration in addition to random seed. Report the discordant fraction as the honest cost of gauge freedom. A rigorous per-cell gauge ensemble — constructed by adding analytically divergence-free perturbations to a fixed drift — would isolate gauge from model-spec variance and is a natural direction for future work; the current library conflates the two.

For cross-condition tests, the depth-matched permutation pipeline is available but should be applied to datasets with a genuine cohort variable (young/old, WT/KO, treatment). Pseudotime terciles do not test the intended non-redundancy hypothesis because both pseudotime and regime derive from the same trained drift field.

## 5. Availability

Software: `github.com/manarai/scJDO` (`scjdo.atlas` module). Install: `pip install scjdo`. Tutorial: `examples/scopatlas/tutorial.ipynb` (synthetic-data validation, ~5 minutes). Real-data application: `examples/scopatlas/paul15_test.ipynb` (Paul15, ~15 minutes CPU). Test suite: `tests/test_scopatlas.py` (48 tests, 2 acceptance gates). Documentation: `docs/scopatlas/`. License: MIT.

---

## Figure

**Figure 1. scOpAtlas on mouse hematopoiesis: replicated non-hyperbolic finding, low-drift-artifact diagnostic, and gauge uncertainty.**

**(A) The replicated $\lambda_{\max}^+ > 0$ observation.** Histograms of per-cell $\lambda_{\max}^+$ on Paul15 (n=2,321, blue) and Setty marrow (n=4,142, orange) under 5-member ensemble primary members trained via `sjd.tl.fit_drift` (hidden=160, depth=3, 2000 epochs). Vertical dashed line at zero marks the hyperbolicity boundary. **Fraction of cells with $\lambda_{\max}^+ < 0$: 0.0% Paul15, 2.9% Setty.** Both distributions bulk above zero. On two canonical mouse hematopoiesis datasets, scJDO-trained drift fields produce essentially no fully-hyperbolic cells. Under the absolute-hyperbolicity classifier, the four-regime partition collapses; under the rank-mode classifier (see panel C), the partition depends entirely on data-driven secondary thresholds. Section 3.2.

**(B) Plasticity indexes low-drift regions.** Per-cell plasticity $P = |\{j: |\lambda_j| < \varepsilon\}|/d$ vs local drift magnitude $\|f(x, t)\|$ (log-scale x-axis). Both datasets show a moderate-to-strong negative Spearman correlation (Paul15 ρ = −0.48, Setty ρ = −0.31). Cells labelled *plastic* sit in low-drift regions of latent space — where the trained network output is close to zero. The label is a low-drift-region indicator, not a distinct spectral state. Users should not interpret *plastic* as a biological plasticity claim without a separate control. Section 3.3.

**(C) Ensemble concordance on Paul15 (rank-mode classifier).** Per-cell concordance distribution across the 5-member ensemble (varying seed × velocity-prior configuration). Dashed line = 0.5 discordance threshold. Mean concordance 0.64; **18.7% discordant at 0.5; 93.4% not fully unanimous** across the ensemble. Honest gauge-uncertainty cost — well below the 33% alarm line, but a real caveat on any single-model regime call. Section 3.4.

---

## Supplementary figures

**Figure S1. Synthetic method validation (3-panel composite).** (Left) Regime classification on a 300-cell 4D synthetic trajectory with hand-designed nonlinear drift produces the intended regime partition. (Middle) 5-member ensemble concordance on the synthetic fixture: mean 0.68, 14% discordant at threshold 0.5. (Right) Random condition assignment on synthetic data produces small non-significant Cramér's V on both all-cells and concordant-only — the acceptance-critical negative control for pipeline calibration on the null.

**Figure S2. Paul15 regime UMAP under rank-mode classifier.** Left: 2,321 Paul15 cells on UMAP colored by rank-classifier regime assignment (four regimes populated: stable, plastic, unstable, deeply-stable). Right: same UMAP colored by DPT pseudotime with root at `7MEP` for spatial reference.

**Figure S3. Paul15 marker cross-validation under rank-mode classifier.** Top row: UMAP overlays of canonical hematopoiesis markers Gata1 (erythroid), Klf1 (terminal erythroid), Elane (myeloid). Bottom row: same UMAP colored by regime. Erythroid TFs concentrate in *plastic* and *deeply-stable* regions; myeloid marker in *stable* region. **Caveat**: with the low-drift-artifact diagnostic from Fig. 1B, the *plastic* label overlaps with terminal-erythroid low-drift regions — the marker concentration is consistent with a low-drift artifact hypothesis as much as with a distinct plasticity state.

**Figure S4. Regime × pseudotime-stage per lineage on Paul15.** Grouped bars with 95% bootstrap CIs. Three lineages significant at $q < 0.05$: erythroid, GMP, myeloid. MEP and Mk non-significant (narrow within-lineage pseudotime range). **Caveat**: pseudotime and regime both derive from the same trained trajectory, so this is a within-model consistency check rather than an orthogonal-cohort test (see §4).

**Figure S5. Paul15 threshold sensitivity (rank-mode classifier).** Regime-fraction heatmaps across a 7×7 grid of $\tau_{\text{upper}}$ (unstable threshold) and $c_{\text{deep}}$ (deep-stable cutoff). Dashed lines mark the auto-resolved defaults. Regime fractions vary smoothly with no sharp transitions near the operating point.

**Figure S6. All cells vs. ensemble-concordant on Paul15.** Per-lineage Cramér's V effect size on all cells (blue) vs. cells with regime concordance ≥ 0.5 (green). Three significant lineages preserve or slightly increase their effect sizes on concordant-only cells — the within-model effect is not primarily a gauge artifact.

---

## References

*(See `references.md` for the bibliography. Verify DOIs before submission. Citations are numbered here for App-Note format.)*

**[1]** Qiu et al. *dynamo*: Mapping transcriptomic vector fields of single cells. Cell 2022.
**[2]** Bocci et al. *SpliceJAC*: transition genes and state-specific gene regulation. Bioinformatics 2022.
**[3]** *scMomentum* — cell fate probabilities from Jacobian analysis (cite specific version).
**[4]** Weinreb et al. Fundamental limits on dynamic inference from single-cell snapshots. PNAS 2018.
**[5]** Wang et al. Multi-scale stochastic dynamics for single-cell transition mapping. Nature Comm 2023.
**[6]** Paul et al. Transcriptional heterogeneity and lineage commitment in myeloid progenitors. Cell 2015.
**[7]** Setty et al. Characterization of cell fate probabilities in single-cell data with Palantir. Nat Biotechnol 2019.

Additional references in `references.md`: La Manno 2018 (RNA velocity), Bergen 2020 (scvelo), Lange 2022 (CellRank), Redd 2026 (scJDO, in prep), Wolf 2018 (scanpy), Benjamini-Hochberg 1995, Miller-Madow 1954, Cramér 1946, Rabin 2012 (sliced Wasserstein).
