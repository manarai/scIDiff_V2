# Figure legends — scJDO manuscript v56

Journal-style legends for every main and supplementary figure, drafted from
the manuscript text and the panel titles embedded in each PDF. Edit for
your target journal's style (word limits, italicisation, panel-letter
formatting).

---

## Figure 1 | scJDO framework overview *(schematic — no code-generated PDF)*

**Overview of the scJDO framework for inferring time-varying dynamical operators from single-cell data.** The pipeline transforms static scRNA-seq measurements into an operator-centric representation of local dynamical sensitivity in four stages. **(a)** Stage 1 — Drift-field estimation: a neural drift field $f_\theta(x, t)$ is learned from a latent representation of cell state via diffusion score matching combined with a residual deterministic component and an optional weak RNA-velocity guidance term. **(b)** Stage 2 — Local operator inference: per-cell Jacobians $J(x, t) = \partial f/\partial x$ are computed by automatic differentiation, summarising how small perturbations are amplified or buffered around each state. **(c)** Stage 3 — Temporal tensor construction: Jacobians are aggregated along pseudotime with Gaussian kernel weighting into a temporal tensor $J \in \mathbb{R}^{T \times d \times d}$. **(d)** Stage 4 — Archetype decomposition: the temporal tensor is unfolded and factorised via semi-nonnegative matrix factorisation (semi-NMF; SVD fallback), yielding a small set of operator archetypes $A_k$ and time-varying activation profiles $a_k(t)$ that reveal when regimes become active, sequentially hand off, or recede.

---

## Figure 2 | Synthetic benchmarks validate the operator machinery and motivate Factor Analysis

**Synthetic benchmarks validate the mathematical core of scJDO, expose a decomposition-choice property on concurrent activations that a shipping opt-in addresses, and a latent-space benchmark motivates Factor Analysis as the default representation.**

*Note*: In the reviewer package, `Figure2.pdf` currently contains the peak-timing null sub-panel (see below). The full manuscript Fig 2 is a composite assembled from the following sub-panel PDFs — recombine in Illustrator:

| Sub-panel source PDF | What it shows |
|---|---|
| `figures/panels_and_diagnostics/fig2_synthetic_benchmark_composite.pdf` | fig2a — 4-panel synthetic benchmark composite (double-well drift-field, dominant-eigenvalue tracking, archetype decomposition, Schrödinger-bridge W₂ test) |
| `figures/main/Figure2.pdf` | fig2b — peak-timing coordination metric validation (24 seeds × 4 synthetic systems) |
| `figures/panels_and_diagnostics/fig2c_embedding_benchmark_8method.pdf` | fig2c — 8-method latent-representation benchmark on Palantir/Setty marrow |
| `figures/panels_and_diagnostics/fig2e_fa_benchmark.pdf` | fig2e — FA-only branch analysis (Ery / DC / Mono top-gene panels) |
| `figures/panels_and_diagnostics/fig2f_fa_pca_scvi_benchmark.pdf` | fig2f — FA vs PCA vs scVI vs LDVAE vs PLS comparison on marrow |

**Suggested panel legend for the composite Fig 2:**

**(a)** *Double-well drift-field recovery.* scJDO learns a smooth continuous-time surrogate of the underlying SDE from noisy trajectory snapshots (endpoint MSE = 0.0005; all simulated trajectories assigned to the correct attractor).

**(b)** *Dominant-eigenvalue tracking.* Along the analytic trajectory the estimated $\mathrm{Re}(\lambda_{\max})$ recovers the ground-truth eigenspectrum (correlation 0.787) and the stable→unstable sign transition.

**(c)** *Archetype-decomposition benchmark.* On a temporal Jacobian sequence with three regimes (stable, transitional, oscillatory), the pipeline recovers three distinct archetype peaks (reconstruction error 0.228).

**(d)** *Schrödinger-bridge benchmark.* Endpoint-constrained transport recovers the analytic target to $W_2 = 0.065$.

**(e)** *Peak-timing coordination metric — synthetic positive control.* Real target metric (◆, 24 seeds) versus 200 null replicates across four synthetic systems (1a sharp handoff, 1b gradual handoff, 2 concurrent overlapping activations, 3 no structure). The metric separates ground-truth-sequential systems (1a, 1b) from concurrent (system 2) and null (system 3) with empirical p < 0.005.

**(f)** *Sequential-recovery Kendall τ across seeds.* Matched Kendall τ on sharp-step (system 1a) = 0.85 with 17/24 seeds at τ ≥ 0.9; on gradual sigmoid (system 1b) = 0.97 with 22/24 seeds at τ ≥ 0.9; per-archetype cosine similarity 0.81 and 0.87 respectively.

**(g)** *Concurrent-activation identifiability.* On overlapping-Gaussian ground truth (system 2), the unregularised semi-NMF returns matched concurrent-peak-pair prevalence 0.025 (vs ground-truth 1.00) because random-init cold starts land in a spread-peaks basin at residual 0.1105, only 0.4 % worse than the GT-init residual 0.1101 (ratio 1.003). An opt-in `tv_lambda` penalty recovers about half the signal (0.49 at $\lambda = 10$) at the cost of destroying sequential recovery on systems 1a/1b; the two operating points are documented in Failure Modes.

**(h)** *Latent-representation benchmark.* Eight embeddings (PCA, FA, ICA, TruncatedSVD, DiffMap, LDVAE, scVI, PLS) compared on the Palantir/Setty marrow dataset under identical pseudotime and branch masks. FA leads the composite score (0.729); DiffMap returns no instability genes and is unscorable on specificity; PLS is supervised and shown as an upper bound.

*Datasets*: Synthetic (analytic ground truth generated in-notebook). Palantir/Setty marrow (4,142 cells; Setty et al. 2019 [1]).
*Source*: `code/figure2_synthetic_and_benchmarks/fig2a_synthetic_benchmark.ipynb` + `fig2b_peak_timing_null_analysis.ipynb` + `fig2c_embedding_benchmark.ipynb` + `fig2e_fa_benchmark.ipynb` + `fig2f_fa_pca_scvi_benchmark.ipynb`; concurrent-activation diagnostics in `concurrent_fix_tests.py`, `init_at_gt_probe.py`, `tv_robustness_test.py`.

---

## Figure 3 | scJDO reveals branch-specific local sensitivity in bone marrow hematopoiesis

**scJDO recovers branch-specific, non-uniform patterns of local sensitivity in the Palantir/Setty CD34+ bone marrow dataset.** Using a 30-dimensional Factor Analysis latent space and Palantir pseudotime, drift fields were fit separately per lineage; each branch captures its trajectory with $R^2 \geq 0.998$.

**(a)** *Palantir trajectory — bone marrow (4,142 cells).* UMAP layout coloured by Palantir pseudotime, with the three terminal lineages annotated: Ery (erythroid), DC (dendritic-cell), Mono (monocyte) and the progenitor compartment.

**(b)** *Local sensitivity — Ery vs DC vs Mono.* $\mathrm{Re}(\lambda_{\max})$ trajectories per branch. Ery (n = 1,151; R² = 0.998) shows a transient early sensitivity regime (peak τ ≈ 0.02, peak Re(λ) ≈ 0.04). DC (n = 1,903; R² = 1.000) shows a near-terminal peak (τ ≈ 0.98). Mono (n = 2,008; R² = 1.000) remains below zero throughout — this branch does not exhibit a positive-eigenvalue sensitivity regime.

**(c)** *Ery archetype activation profiles (K = 5).* Peak order A1 → A4 → A5 → A3 → A2. The archetype anchored at τ ≈ 0.02 captures the transient erythroid maturation program.

**(d)** *DC archetype activation profiles (K = 5).* Peak order A1 → A4 → A5 → A3 → A2, with the late-terminal archetype (A2 at τ ≈ 0.98) coinciding with the sensitivity peak.

**(e)** *Top instability-driving genes per branch.* Ery top: HBB, AHSP, CA1, FAM178B, APOC1 — canonical erythroid maturation markers; GATA1 recovered without supervision as the top transcription-factor regulator. DC top: TK1, JCHAIN, IRF8, C12ORF75, CLSPN — a mixed immune/DC-identity + proliferative signature. Mono returned no positive-eigenvalue instability genes under the same criterion, consistent with a comparatively buffered operator trajectory in this lineage.

**(f)** *Mono archetype activation profiles (K = 5)*, shown for reference; peak order A1 → A4 → A3 → A2 → A5.

*Dataset*: Palantir/Setty CD34+ bone marrow (4,142 cells; Setty et al. 2019 [1]).
*Source*: `code/figure3_hematopoiesis/figure3_hematopoiesis.ipynb`. Palantir-vs-DPT ordering agreement quantified in `palantir_dpt_agreement.py` + `palantir_dpt_agreement_extra_variants.py` (branch-restricted Spearman ρ range 0.839–0.988 across nine tested configurations).

---

## Figure 4 | scJDO is concordant with Dynamo but distinct from splicing-derived Jacobians

**scJDO's operator-level structure is concordant with Dynamo at the gene and eigenvector levels but distinct from SpliceJAC on the Setty 2019 CD34+ bone marrow dataset.** The concordance pipeline operates on the 6,482 filtered genes shared by all three methods and uses top-30 transition genes per cluster (10 clusters).

**(a)** *Cluster-level max Re(λ) — scJDO vs SpliceJAC / Dynamo.* Per-cluster maximum real eigenvalue for each method, scattered against scJDO. Systematic ranking inversion vs Dynamo (Spearman ρ = −0.818, p = 3.8 × 10⁻³) and vs SpliceJAC (ρ = −0.77) — indicated by the annotated ρ values. scJDO's peak-instability cluster is CLP; Dynamo's is Precursors.

**(b)** *Leading unstable-eigenvector direction agreement.* Per-cluster |cos| between the leading unstable eigenvector from scJDO and from each comparator. Dynamo mean |cos| = 0.100 (above null 95th percentile in 8 of 10 clusters, 80 %); SpliceJAC mean |cos| = 0.030 (above null in 5 of 10, 50 %). The empirical null 95th percentile is ≈ 0.024.

**(c)** *Transition-gene overlap — mean Jaccard at top-30.* Dynamo mean Jaccard 0.0749; SpliceJAC 0.0051. Under a random-gene null over the shared 6,482-gene set, Dynamo overlap is 32× the null and SpliceJAC 2.2×. Under the more conservative null restricted to the shared top-2,000 HVG vocabulary that both methods actually draw from, Dynamo overlap is ~10× and SpliceJAC ~0.7× — i.e., SpliceJAC's overlap is at or below the analytic null and not distinguishable from random draws over the shared HVG set.

*Dataset*: Setty 2019 CD34+ bone marrow, 5,780 cells, 6,482 shared filtered genes; auto-downloaded via `scvelo.datasets.bonemarrow()`.
*Source*: `code/figure4_concordance/fig4_04_concordance_metrics.ipynb` (primary pipeline; produces `concordance_metrics.json`); `top_hvg_null.py` computes the two-null fold decomposition from the JSON. Per-method operator computations in `fig4_01_scjdo_operators.ipynb`, `fig4_02_splicejac_operators.ipynb`, `fig4_03_dynamo_operators.ipynb`.

---

## Figure 5 | scJDO reveals distinct operator handoffs during MEF-to-iPSC reprogramming in SCP295

**scJDO separates productive transition-state instability in the MET corridor from a distinct diverted stromal instability program in a densely sampled iPSC reprogramming time course.** Because SCP295 supplies an *observed* temporal axis (experimental day) rather than an inferred pseudotime, this is the setting in which the manuscript makes timing-resolved claims.

**(a)** *Reprogramming trajectory — SCP295 (24,999 cells, days 0–18).* Force-directed layout coloured by experimental day, with the three fate-annotated compartments (IPS, Stromal, MET) and the MEF/Progenitor starting state.

**(b)** *Local sensitivity — IPS / Stromal / MET / IPS_late.* $\mathrm{Re}(\lambda_{\max})$ trajectories per branch. MET (n = 5,267; R² = 0.996) shows a productive-instability peak at τ ≈ 0.77 within the fitted branch model. Stromal (n = 6,350; R² = 1.000) shows a distinct late sensitivity regime (τ ≈ 0.98). IPS (n = 5,249; R² = 1.000) shows a modest early sensitivity peak (τ ≈ 0.02). IPS_late (n = 8,454; R² = 1.000) shows a weaker terminal proliferative regime (τ ≈ 0.98). Because absolute eigenvalue magnitudes are not comparable across independently fitted branch drift fields, we do not rank instability magnitudes across branches.

**(c–f)** *Per-branch archetype activation profiles (K = 5).* Productive branches (IPS, MET) execute a sequential handoff from an early MEF-exit archetype (A1, peak τ ≈ 0.02 on the MET corridor) through intermediate archetypes to a late pluripotency/metabolic-shift archetype (A2 on MET at τ ≈ 0.61; A3 on IPS at τ ≈ 0.69 under TV = 0.1). The Stromal-diversion branch executes no such completion: it maintains the early MEF-exit archetype and instead activates a distinct stress-associated archetype.

**Top instability genes and regulators per branch** (from the per-branch panels in `panels_and_diagnostics/fig5_instab_*.pdf`):
- MET: *Tagln, Ccl7, Id3, Gpx3, Acta2*; regulators *Trp53, Sp1, Nfkb1, Smad3, Rela*.
- Stromal: *Cdkn2a, Gas6, Lmo4, Mxra8*; NFκB/p53-associated regulators.
- IPS (early): *Col6a3, Col3a1, Col5a1, Col1a1, Col6a1* — fibroblast/ECM exit.
- IPS_late: *Birc5, Cenpf, Top2a, Ube2c, Cdca3* — cell-cycle expansion of successfully reprogrammed cells.

**Companion `Figure5_TV01.pdf`** — the same pipeline re-run with a mild total-variation prior on semi-NMF activations (`TV_LAMBDA = 0.1`). On the MET branch the archetype peak order is byte-identical to TV = 0 (A1 (0.02) → A4 → A2 (0.61) → A5 (0.67) → A3 (≈ 0.75)); the A1→A2 crossover at t ≈ 0.69 in productive branches is therefore not a decomposition artifact of the unregularised default.

*Dataset*: SCP295 reprogramming (Schiebinger et al. 2019 [2]) — 24,999 cells, days 0–18, obs['day'] + obs['cell_set'].
*Source*: `code/figure5_reprogramming/figure5_reprogramming.ipynb` (TV = 0) + `figure5_reprogramming_tv01.ipynb` (TV = 0.1). Bandwidth sensitivity (peak τ 0.60 → 0.95 as h sweeps 0.01 → 0.15) documented in `code/supp_bandwidth_sweep/bandwidth_sweep_dc_met.py`.

---

## Figure 6 | Per-perturbation Schrödinger bridges recover coherent regulatory programs in K562 CRISPRi Perturb-seq

**scJDO's Schrödinger-bridge formulation extends the framework to perturbational data where no natural pseudotime exists.** Forward bridges from non-targeting control (NT) cells to each of three lncRNA-knockdown populations (PVT1, MALAT1, PSMA3-AS1) recover perturbation-specific gene programs; MALAT1 and PSMA3-AS1 converge on a shared, multiple-testing-robust G2M-checkpoint arrest program.

**(a)** *K562 Perturb-seq — NT + 3 lncRNA targets.* UMAP layout with NT (n = 1,238; grey) and the three targets (PVT1 n = 135; MALAT1 n = 126; PSMA3-AS1 n = 89) coloured by perturbation identity.

**(b)** *Forward-bridge sensitivity — per perturbation.* $\mathrm{Re}(\lambda_{\max})$ along bridge time (NT → perturbed) for each target. All three peak near the perturbed endpoint (peak_t ≈ 0.97). PVT1 max λ = +3.170 (strongest raw magnitude but see validation panel); MALAT1 = +0.311; PSMA3-AS1 = +0.644.

**(c–e)** *Per-perturbation archetype activation profiles (K = 5).* PVT1, MALAT1, PSMA3-AS1 respectively. Order and peak τ per archetype annotated on each panel.

**Top forward-instability gene and regulators per perturbation:**
- PVT1: top gene *HNRNPA2B1*; top regulator *SP1*. Reported as a gene-level nomination only — pathway enrichment (PI3K–AKT–mTOR / mTORC1) rested on a single overlapping gene and did not survive multiple-testing correction (FDR ≈ 0.64).
- MALAT1: top gene *TOP2A*; top regulator *TP53*.
- PSMA3-AS1: top gene *TOP2A*; top regulator *TP53*.

*Dataset*: K562 CRISPRi Perturb-seq (10x Genomics `84K_K562_5pv3_GEMX_CRISPR_aggregate`, CC BY 4.0). Replogle 2022 [13] genome-scale lncRNA sublibrary.
*Source*: `code/figure6_k562_perturbseq/figure6_k562.ipynb`.

## Figure 6 — Validation | Hallmark enrichment + permutation control

**Statistical support for the Fig 6 gene-program claims.**

**(a)** *Hallmark over-representation in scJDO top-50 genes (−log₁₀ FDR).* MALAT1 G2M-Checkpoint FDR = 2 × 10⁻¹⁷ with 12 overlapping mitotic genes (AURKA, CCNA2, CDC20, CENPE, CENPF, HMMR, KIF2C, MKI67, PLK1, TOP2A, TPX2, …). PSMA3-AS1 G2M-Checkpoint FDR = 8 × 10⁻¹³ with 10 overlapping mitotic genes. PVT1 shows no G2M enrichment; its top hallmark pathways (MTORC1, PI3K-AKT-MTOR) rest on single overlaps that do not survive FDR correction.

**(b)** *Permutation control (PVT1, n = 20).* Null distribution of the forward-bridge $\mathrm{Re}(\lambda_{\max})$ under permuted perturbation labels; observed λ = +3.17 corresponds to empirical p = 0.143 — not distinguishable from random. This permutation was run for PVT1 only (the strongest raw signal); MALAT1 and PSMA3-AS1 significance rests on the pathway-level enrichment in panel (a), not on the eigenvalue magnitude. Consistent with our caveat that operator magnitudes are uninformative in this static perturbation regime — the value of the analysis is in the recovered program content, not in the instability magnitude.

*Source*: `code/figure6_k562_perturbseq/figure6_k562.ipynb` (same notebook as Fig 6).

---

## Supplementary Figures

### Fig. S1 | Multiome joint RNA + ATAC FA integration

**Joint FactorAnalysis integration of a 10x Multiome (RNA + ATAC) dataset for the Supplementary Note.** Shared 30-D FA latent space fit on 2,000 top-variable RNA HVGs and 10,000 top-variable ATAC peaks (post-TF-IDF), with per-factor RNA vs ATAC vs joint contribution shown per component. Integration QC (per-cell neighbourhood overlap between RNA-only, ATAC-only, and joint latents) is in Fig. S2.

*Source*: `code/supp_multiome_note/supp_multiome_FA_integration.ipynb`.

### Fig. S2 | Multiome integration QC

**Per-cell neighbourhood-overlap and clustering-similarity metrics** confirming that the joint FA latent space preserves both modality-specific neighbourhood structures.

*Source*: same notebook as Fig. S1.

### Fig. S3 | Multiome marker-ordering diff before vs after cell-cycle regression

**Marker-gene pseudotime-ordering shift before and after regressing S/G2M cell-cycle scores from the latent space.** In proliferative tissue, cell-cycle-associated variance dominates the raw drift field; regression aligns the ordering to lineage biology.

*Source*: `code/supp_multiome_note/supp_multiome_drift_cellcycle_regressed.ipynb`.

### Fig. S4 | Multiome ExcNeuron instability — cell-cycle-regressed

**Instability trajectory for the excitatory-neuron branch after cell-cycle regression.**

*Source*: same notebook as Fig. S3.

### Fig. S5 | Multiome ExcNeuron instability — uncorrected

**Same excitatory-neuron branch analysed without cell-cycle regression**, for contrast. The uncorrected version overweights cell-cycle-associated genes (Cdc20, Ccnb1, etc.) at the expense of the lineage program.

*Source*: `code/supp_multiome_note/run_uncorrected_multiome_pipeline.py`.

### Fig. S6 | Bandwidth sweep — DC / MET

**Supported peak pseudotime as a function of the temporal-kernel bandwidth $h$**, on a simplified reproduction of the MET-branch pipeline (SCP295, 3,344 MET-subset cells, 1,000 training epochs, no velocity-guidance prior). The supported peak τ shifts monotonically from 0.60 at h = 0.01 to 0.95 at h = 0.15; only h = 0.01 returns a positive peak eigenvalue on this simplified pipeline. Included so that any downstream timing claim can be audited against a bandwidth-sweep sensitivity check on the same pipeline.

*Source*: `code/supp_bandwidth_sweep/bandwidth_sweep_dc_met.py` + `make_supp_fig_bandwidth.py`.

### Fig. S7 | Saddle-localization label-blind test

**Test whether the flow topology of the trained drift field encodes the analytic saddle location on the printed toggle-switch bifurcation system.** Two label-blind signatures — drift-magnitude collapse and forward-trajectory divergence — computed from the same trained drift field that produced the failed $\lambda_{\max}(\tau)$ recovery, and scored against a pre-committed window around the corrected analytic pitchfork ($\tau_{\mathrm{crit}} \approx 0.004$, window $[0, 0.054]$). Neither signature localises the saddle, supporting the interpretation that the density-dominance limit is upstream of the scalar readout rather than a peculiarity of $\mathrm{Re}(\lambda_{\max})$.

*Source*: `code/supp_saddle_localization/saddle_localization_v3_label_blind.ipynb`. Companion peak-value analyses in `saddle_readouts.py` and `saddle_dense_sampling.py`.

---

## Illustrator workflow notes

- **Fonts**: matplotlib saves fonts as embedded TrueType (`pdf.fonttype = 42`) so text stays selectable/editable in Illustrator.
- **Bundle for opening**: `figures/all_figures.pdf` (14-page concatenation of all main + supplementary PDFs) can be opened in Illustrator's *Open* dialog, which will prompt for a page number.
- **Individual editing**: open the per-file PDFs directly from `figures/main/` or `figures/supplementary/`.
- **After editing**: re-save the edited PDF back to its original location; regenerate the matching PNG at 300 DPI with

    pdftoppm -r 300 -png -singlefile figures/main/Figure3.pdf figures/main/Figure3
