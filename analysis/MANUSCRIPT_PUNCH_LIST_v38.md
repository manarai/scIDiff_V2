# Manuscript v38 punch list — what to change so the paper matches the code

All items below were verified against the actual notebooks, saved CSVs/JSONs, and re-runs of the missing analyses.
The **Trust the CODE** direction was applied — every recommendation below fixes the manuscript, not the pipeline.

Fixes already applied to notebooks by this pass are noted at the bottom.

---

## 1 — Figure 2 (synthetic benchmarks + FA embedding)

| # | Manuscript v38 says | Code / notebook actually shows | Recommended manuscript change |
|---|---|---|---|
| 1.1 | Sigmoid handoff (System 1b): matched Kendall τ = 0.97, 22/24 seeds ≥ 0.9; sharp-step (1a): τ = 0.85, 17/24 seeds ≥ 0.9 | CSV matches manuscript exactly (`decomposition_recovery_metrics.csv`: 0.9667, 22/24; 0.8529, 17/24). | Keep — no change. |
| 1.2 | Concurrent (System 2): residual 0.1101 (GT concurrent) vs 0.1105 (spread peaks), ratio 1.003 | New `init_at_gt_probe.json` (8 seeds): warm-start-at-GT ALS mean residual 0.1101; cold-start ALS mean 0.1105 (see `analysis/results/figure2/init_at_gt_probe.json`). Concurrent-peak-pair prevalence = 1.0 on all 8 seeds under warm start. | Keep — no change. Ratio 1.003 confirmed. Static W=A_gt residual (0.1103) is a separate number and shouldn't be conflated with the ALS residuals. |
| 1.3 | tv_lambda = 10 → 0.49 concurrent prevalence | Confirmed in existing CSV. | Keep. |
| 1.4 | Composite scoring (Supplementary Table 1) — with timing: PLS 0.622; without timing: FA 0.729, PLS 0.732, FA-second 0.706 | Confirmed in `benchmark_scores.csv` + `fa_composite_reanalysis.csv`. | Keep. |
| 1.5 | Ery ∩ DC top-N Jaccard = 0 | Confirmed — but the shipped CSV stores top-15, notebook stdout says "top-30" in one place and the manuscript quotes "top-20". | Pick one N and use it consistently in text, notebook stdout, and the saved CSV filename. Top-20 (the manuscript's choice) works — expand the saved CSV. |

## 2 — Figure 3 (Paul15 hematopoiesis)

| # | Manuscript v38 says | Code / notebook actually shows | Recommended manuscript change |
|---|---|---|---|
| 2.0 | **"The Paul15 hematopoiesis dataset is available via the Scanpy package."** (Data Availability) and "Paul15 hematopoiesis" in Results text | The dataset actually loaded by `Figure3_FA.ipynb` cell 4 is **not** `sc.datasets.paul15()`. It is the **4,142-cell mouse bone-marrow sample from the Palantir tutorial** (`examples/marrow_sample_scseq_counts.h5ad`), originally from Setty *et al.* 2019 CD34+ profiling. The n=1151/1903/2008 branch counts derive from this dataset, not from Paul et al. 2015. | **Rename throughout the manuscript.** Suggested wording — Data Availability: "The mouse bone-marrow sample dataset used for Figure 3 is the 4,142-cell dataset from the Palantir tutorial (from Setty *et al.* 2019); available at `https://github.com/dpeerlab/Palantir` under `data/marrow_sample_scseq_counts.h5ad`." Results text: replace "Paul15 hematopoiesis" with "the Palantir/Setty mouse bone-marrow sample" (or shorter: "the Palantir mouse-marrow sample"). Also cite Setty *et al.* Nat Biotech 2019 (already ref 2) rather than Paul et al. 2015 (ref 1) for this dataset. |
| 2.1 | Ery n=1151, DC n=1903, Mono n=2008; R² = 0.998/1.000/1.000 | Confirmed in `branch_summary.csv`. | Keep (with corrected dataset name from 2.0 above). |
| 2.2 | Ery top instability genes include HBB, AHSP, CA1, HBD; top regulator GATA1 | Confirmed in `instability_genes_Ery.csv` and `regulators_Ery.csv`. | Keep. |
| 2.3 | DC signature JCHAIN, IRF8, IRF7, PLD4, RNASE6, TK1, RRM2 | Confirmed — all seven in top-10 of `instability_genes_DC.csv`. | Keep. |
| 2.4 | Mono: no positive-eigenvalue regime, no instability genes / regulators | Confirmed — `instability_genes_Mono.csv` empty, `regulators_Mono.csv` header-only, `max_eig = −0.0417`. | Keep. |
| 2.5 | **"Palantir-vs-DPT branch-restricted rank correlation ρ = 0.08 on the erythroid branch"** (used to justify declining a peak-timing claim) | ✗ **DOES NOT REPRODUCE under any tested configuration.** Nine variants across two scripts, ρ range **[0.839, 0.988]** — nothing anywhere near 0.08: <br>• Global Palantir & DPT restricted to Ery, k=15 FA:                                 ρ = **0.842** <br>• DPT rerun on Ery subset (k=15 FA), matched local root:                              ρ = **0.972** <br>• Both Palantir & DPT rerun on Ery subset:                                            ρ = **0.961** <br>• DPT on PCA (k=15, global root):                                                     ρ = **0.874** <br>• DPT on PCA (k=15) restricted to Ery, local root:                                    ρ = **0.839** <br>• DPT on FA with scanpy-default k=30, global root:                                    ρ = **0.841** <br>• DPT on FA with 5 randomly-sampled early-progenitor roots (Palantir pt < 0.05):     ρ ∈ **[0.966, 0.988]** <br>• DPT on raw HVG log expression (no dim-reduction) — see `ery_palantir_dpt_extra_variants.json` for D6 result <br>DC and Mono all likewise land ρ ≥ 0.83 under every variant. The full JSON dumps are at `analysis/results/figure3/ery_palantir_dpt_agreement.json` and `.../ery_palantir_dpt_extra_variants.json`. | **Retract the ρ = 0.08 claim.** After 9 variants spanning FA vs PCA vs raw expression, k=15 vs 30, global root vs local root vs random-early-cell root, nothing lands near 0.08. Either the number came from a bug in an earlier iteration or from a configuration I haven't been able to reconstruct. Recommended manuscript text: replace with the actual measured value (Palantir vs DPT on Ery, k=15 in FA, matched root: ρ = 0.842) and rewrite the "we withhold peak-timing on Ery because Palantir vs DPT disagree" reasoning — that reasoning is not supported by the actual data. |
| 2.6 | "absolute λ_max thresholds classifying modes: stable ≤ −0.05, sensitive ≥ +0.05, otherwise plastic" | Notebook code only implements the `+0.05` sensitive cutoff (`n_sens = (eig > 0.05).sum()`) — the `−0.05` stable cutoff is not tested anywhere. With `max_eig` = 0.0438 (Ery) and 0.0420 (DC), the `+0.05` rule reports `n_sensitive = 0/200` for BOTH branches — inconsistent with the manuscript's characterization of Ery and DC as the sensitive branches. | Either (a) lower the sensitive cutoff to +0.04 so the classification matches the biological reading of Ery/DC as sensitive; or (b) drop the specific ±0.05 numbers from the Methods paragraph and describe the classification qualitatively. The +0.05 threshold as currently stated puts Ery and DC in the "plastic" bin, which contradicts the Results narrative. |

## 3 — Figure 4 (Setty CD34+ concordance)

| # | Manuscript v38 says | Code shows | Change |
|---|---|---|---|
| 3.1 | Fold-over-null: Dynamo 161× (random), 49× (top-HVG); SpliceJAC 29.8× / 9.1× | `concordance_metrics.json` = 160.91, `top_hvg_null.json` = 49.39; SpliceJAC 29.77 / 9.14 | Keep. |
| 3.2 | Mean \|cos\| leading unstable eigenvector: Dynamo 0.23 (100% > null 95th), SpliceJAC 0.05 (40%) | Confirmed in JSON. | Keep. |
| 3.3 | Branch-level ρ = −0.903, p = 3.4 × 10⁻⁴; scJDO argmax = CLP, Dynamo = Mega | Confirmed. | Keep. |
| 3.4 | "both methods draw from the same top-2000 HVG vocabulary" | The pipeline (`fig4_01/02/03.ipynb`) works over the 6,482 filtered genes and doesn't force a 2000-HVG intersection. The top-2000 HVG is only used inside `top_hvg_null.py` as the null's universe. | Rephrase: "the conservative null restricts the random draw to the shared top-2000 HVG universe that overlaps the working gene set" — i.e., it's a null construction, not a pipeline gate. |

## 4 — Figure 5 (SCP295 reprogramming)

| # | Manuscript v38 says | Code shows | Change |
|---|---|---|---|
| 4.1 | **n = 27,538 cells** | Notebook `Figure4.ipynb` cell 0 explicitly says "**24,999-cell stratified subsample**"; panel-a title, `branch_summary.csv`, and cell 3 stdout all print n = 24,999. `build_scp295_h5ad.py` (the pipeline that produces the loaded h5ad) generates exactly 24,999 cells. | **Change manuscript to n = 24,999** (author-confirmed as authoritative). Also add the stratified-subsampling detail to Methods: "n = 24,999 cells drawn as a stratified subsample of the SCP295 upload, produced by `Figures_notebook/build_scp295_h5ad.py`." |
| 4.2 | MET top genes Tagln, Ccl7, Id3, Gpx3, Acta2 | Confirmed exactly. | Keep. |
| 4.3 | MET top regulators Trp53, Nfkb1, Smad3, Rela | Actual `regulators_MET.csv` order is **Trp53, Sp1, Nfkb1, Smad3, Rela** (Sp1 is rank 2, silently omitted from manuscript's list). | Add **Sp1** (rank 2) to the MET top-regulator sentence, or note explicitly that Sp1 was excluded on biological grounds. |
| 4.4 | MET bandwidth sweep: peak τ monotonic 0.60 → 0.95 across h=0.01→0.15, only h=0.01 positive | Confirmed exactly in `bandwidth_sweep_dc_met.json` (0.604 → 0.946; only h=0.01 has positive `peak_lam`). Full-pipeline peak τ ≈ 0.77 also confirmed. | Keep. |
| 4.5 | Stromal senescence signature (Cdkn2a, Gas6, Lmo4, Mxra8), IPS Col-family, IPS_late Birc5/Cenpf/Top2a/Ube2c/Cdca3 | All confirmed exactly. | Keep. |
| 4.6 | **A1→A2 crossover at experimental-day-normalized time ≈ 0.69 on MET and IPS** | On the MET branch, the actual archetype peak order is **A1 → A4 → A2 → A5 → A3** (A1 at t≈0.02, A4 at 0.25, A2 at 0.61, A5 at 0.67, A3 at 0.77). The "A1→A2 handoff at 0.69" phrasing doesn't cleanly describe this sequence. On the IPS branch, TV=0 late-archetype (A5) peaks at t≈0.72; under TV=0.1 the late archetype (A3) peaks at t=0.69. | Rephrase to describe the actual archetype sequence. Suggested wording: "productive branches execute a sequential handoff from an early MEF-exit archetype (A1, t≈0.02 on MET) to a late pluripotency/metabolic-shift archetype (A2 on MET at t≈0.61; A3 on IPS at t≈0.69 under TV=0.1)". Do not claim an A1→A2 handoff on IPS — the manuscript is conflating archetype labels across branches. |
| 4.7 | **"TV=0.1 preserves archetype-column order byte-identically"** | Byte-identical column ordering holds on the MET branch only. IPS, Stromal, and IPS_late orderings differ between TV=0 and TV=0.1. | Restrict scope to "MET-branch archetype-column ordering is byte-identical between TV=0 and TV=0.1 (A1→A4→A2→A5→A3 in both); other branches show semi-NMF's expected column-permutation ambiguity under a non-zero smoothness prior, without disturbing the overall temporal structure". |
| 4.8 | Synthetic controls (Systems 1a/1b): matched Kendall τ = 1.00 on 8/8 seeds at λ=0.1 | ✓ **Confirmed exactly.** `analysis/results/figure2/tv_robustness_test.json`: 1a τ_mean = 1.000, 8/8 seeds ≥ 0.9; 1b τ_mean = 1.000, 8/8 seeds ≥ 0.9. The full TV grid also shows why λ=0.1 is the safe operating point: at λ=1.0 System 1b collapses to τ_mean 0.590 (1/8 seeds); at λ=10 both systems fall apart. | Keep. |

## 5 — Figure 6 (K562 CRISPRi Perturb-seq)

| # | Manuscript v38 says | Code shows | Change |
|---|---|---|---|
| 5.1 | MALAT1 G2M FDR ≈ 2×10⁻¹⁷ over 12 mitotic genes; PSMA3-AS1 G2M FDR < 10⁻¹² | Confirmed exactly (MALAT1 = 2.06e-17 / k=12; PSMA3-AS1 = 7.59e-13 / k=10). | Keep. |
| 5.2 | CDC20, CENPE, CENPF, TOP2A in shared MALAT1/PSMA3-AS1 program | Confirmed in `top_genes_per_target.json` and G2M/MITOTIC_SPINDLE overlap lists. | Keep. |
| 5.3 | PVT1 top instability gene = HNRNPA2B1; PVT1 PI3K-AKT-mTOR / mTORC1 rests on 1 gene, FDR ≈ 0.64 | Confirmed exactly (`figure5_summary.csv`, `hallmark_ora.csv`: PI3K_AKT_MTOR k=1 FDR=0.638; MTORC1 k=1 FDR=0.655; both HSP90B1). | Keep. |
| 5.4 | Bridge eigenvalue null: empirical p = 0.143 | Confirmed in `permutation_pvt1.json` (p=0.1428…). Manuscript reads as if this applies to all three perturbations — the null was actually run on **PVT1 only**. | Clarify: "For PVT1 — the strongest signal — the bridge eigenvalue was not distinguishable from a random null (empirical p = 0.143, 20 permutations). We did not run the permutation control on MALAT1 or PSMA3-AS1; their significance rests on the pathway-level enrichment, not on the eigenvalue magnitude." |

## 6 — Supplementary: bandwidth sweep (Figure 5 supplement)

| # | Manuscript v38 says | Code shows | Change |
|---|---|---|---|
| 6.1 | 3,344 MET cells, 1,000 epochs, no velocity prior, 0.60→0.95 monotonic peak-τ shift | Confirmed exactly (`bandwidth_sweep_dc_met.json` + log). | Keep. |

## 7 — Supplementary: multiome cell-cycle Note

| # | Manuscript v38 says | Code shows | Change |
|---|---|---|---|
| 7.1 | E18 mouse-brain multiome, 4,878 cells | Confirmed. | Keep. |
| 7.2 | Post-correction top-9 ExcNeuron = Raly, Gnai2, Ndufa11 + others, no mitotic genes | Confirmed exactly in `instability_genes_cc_ExcNeuron.csv`. | Keep. |
| 7.3 | Canonical developmental marker ordering invariant to correction | Confirmed (all 8 markers unchanged mean pseudotime). | Keep. |
| 7.4 | **"pre-correction top-9 was dominated by Cdc20 / Ccnb1 / Cenpf"** | ✗ **NOT REPRODUCIBLE from current artifacts.** The "uncorrected" input (`adata_multiome_fa.h5ad`) already has all three FA latents (RNA, ATAC, joint) orthogonal to S/G2M scores at the leading components — pre-integration correction has been applied upstream. As a result the actual saved pre-correction top-9 is Rhbdl3, Eomes, Neurog2, Celsr1, Gm29260, Smoc1, Chd7, Mgp, C030005K06Rik — a clean neurogenic signature with zero mitotic hits. Details: `analysis/results/supp_multiome/precorrection_reproducibility_check.json`. | Two options: (a) re-run the multiome pipeline from raw 10x E18 data with per-modality CC regression turned OFF, and persist that output — if it reproduces mitotic dominance, keep the current wording; (b) rewrite the Note to describe what the current pipeline actually shows (the pipeline is robust to further latent-space CC regression; both uncorrected and CC-regressed pipelines return non-mitotic developmental top-9 lists). The current wording is not backed by any file this evaluator could locate. |

## 8 — Supplementary: bifurcation-saddle localization

| # | Manuscript v38 says | Code shows | Change |
|---|---|---|---|
| 8.1 | τ_crit = 0.167; ±0.05 pre-committed window; readouts Re(λ_max), tr(J), P⊥JP⊥, FTLE | Confirmed. | Keep. |
| 8.2 | Baseline per-readout peak τ = 0.83 / 0.95 / 0.05 / 0.85 / 0.78 across **"2 seeds"** | Values match, but `saddle_readouts.py` uses a **single SEED=42**, not 2 seeds. The dense-sampling control uses 2 seeds; the baseline does not. | Change "mean across 2 seeds" to "single seed" for the baseline (or add a second baseline seed to the script). |
| 8.3 | Dense-sampling control: 5× density, 2 seeds, peaks 0.93 / 0.95 / 0.05 / 0.86 / 0.86; 0/2 seeds pass | Confirmed exactly (`saddle_dense_sampling.json`). | Keep. |
| 8.4 | FTLE parameters: horizon 0.5, **100 Euler steps**, δ = 10⁻³ | Both `saddle_readouts.py` and `saddle_dense_sampling.py` use **`steps=50`**, not 100. δ and horizon are correct. | Change "100 Euler steps" → "50 Euler steps" (or update the scripts to `steps=100` if that is the intended value; the qualitative conclusion is unlikely to move). |

---

## Fixes already applied to notebooks by this pass

1. `Figures_notebook/Figure2_synthetic_peak_timing_null_analysis.ipynb` and `analysis/figure2_synthetic_and_benchmarks/fig2b_peak_timing_null_analysis.ipynb` — cell 16 markdown and cell 18 figure suptitle updated from stale (τ=0.81/1.00, 20/24/24/24, cosine 0.86–0.95) to CSV-consistent values (τ=0.85/0.97, 17/24/22/24, cosine 0.81/0.87).
2. `Figures_notebook/Figure5.ipynb` and `analysis/figure6_k562_perturbseq/figure6_k562.ipynb` — cell 19 summary updated: PVT1 "MTORC1 FDR 2×10⁻⁶" (stale) replaced with actual "PI3K_AKT_MTOR single-gene overlap HSP90B1, FDR ≈ 0.64 — does not survive multiple-testing"; MALAT1 "G2M FDR 4×10⁻²⁹, 22-gene overlap" (stale, ran on top-50) replaced with actual "G2M FDR ≈ 2×10⁻¹⁷, 12-gene overlap" (matches manuscript and current CSV); PSMA3-AS1 "G2M FDR 5×10⁻³¹, 23-gene overlap" → "G2M FDR ≈ 8×10⁻¹³, 10-gene overlap"; permutation-control scope clarified as PVT1-only.
3. `analysis/figure2_synthetic_and_benchmarks/init_at_gt_probe.py` — added JSON output (`analysis/results/figure2/init_at_gt_probe.json`) with per-seed residuals + concurrent-peak prevalence over all 8 seeds. Confirms manuscript's 0.1101 vs 0.1105 numbers exactly.
4. `analysis/figure2_synthetic_and_benchmarks/tv_robustness_test.py` — added JSON output (`analysis/results/figure2/tv_robustness_test.json`) with per-seed matched Kendall τ across TV grid {0, 0.01, 0.1, 1.0, 10.0} for Systems 1a and 1b.
5. `analysis/figure3_hematopoiesis/palantir_dpt_agreement.py` (new) — computes three-variant Palantir-vs-DPT branch-restricted Spearman ρ on Paul15 with matched root and k=15 kNN. Output: `analysis/results/figure3/ery_palantir_dpt_agreement.json`.
6. `analysis/results/supp_multiome/precorrection_reproducibility_check.json` (new) — documents the multiome "pre-correction Cdc20/Ccnb1/Cenpf" claim's non-reproducibility from currently shipped artifacts.

## Items requiring human decision (not autonomous edits)

- **Fig 5 cell count 24,999 vs 27,538:** the notebook is authoritative here (24,999). If 27,538 was the correct dataset size, the SCP295 loading step differs from what the notebook does. Decide whether to update the manuscript number or reload data.
- **Fig 3 Palantir-vs-DPT ρ = 0.08:** no variant reproduces this value. Either provide the exact runnable configuration that produced it, or drop the claim.
- **Multiome pre-correction mitotic dominance:** would need a raw-data reprocessing to demonstrate; requires access to the raw ATAC modality which isn't in the h5ad.
