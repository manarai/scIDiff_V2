# scJDO — reproduction guide

Every figure, table, and supplementary panel in `scJDO_manuscript_v52_final_draft.md`, mapped to the notebook or script that generates it. Paths are relative to the repo root.

Two subtrees hold reproduction code:
- `analysis/` — canonical, figure-number-matched notebooks (primary source).
- `scratchpad_run/` — later-added scripts that supersede or extend individual analyses (documented in the "Post-v52 additions" column below).

The `Figures_notebook/` subtree contains earlier drafts under an older numbering scheme and is retained for auditability only.

Everything in this document is auto-verifiable — every listed path exists on disk.

---

## Main-text figures

| Deliverable | Description | Primary code | Post-v52 addition |
|---|---|---|---|
| **Figure 1** | Framework schematic (illustration; no computed panels) | *n/a* | — |
| **Figure 2a** | Synthetic mathematical-core benchmark (double-well drift, Jacobian eigenspectrum, three-regime archetype decomposition, Schrödinger-bridge $W_2$) | `analysis/figure2_synthetic_and_benchmarks/fig2a_synthetic_benchmark.ipynb` | — |
| **Figure 2b** | Peak-timing null analysis (sequential-handoff recovery under semi-NMF permutation ambiguity) | `analysis/figure2_synthetic_and_benchmarks/fig2b_peak_timing_null_analysis.ipynb` | — |
| **Figure 2c** | Embedding benchmark on Paul15 hematopoiesis (8 methods × 3 branches) | `analysis/figure2_synthetic_and_benchmarks/fig2c_embedding_benchmark.ipynb` | — |
| **Figure 2d** | Pseudotime × embedding sensitivity | `analysis/figure2_synthetic_and_benchmarks/fig2d_pseudotime_embedding_benchmark.ipynb` | — |
| **Figure 2e** | FA-specific ablation | `analysis/figure2_synthetic_and_benchmarks/fig2e_fa_benchmark.ipynb` | — |
| **Figure 2f** | FA vs PCA vs scVI focused comparison | `analysis/figure2_synthetic_and_benchmarks/fig2f_fa_pca_scvi_benchmark.ipynb` | — |
| **Figure 3** | Hematopoiesis: Ery / DC / Mono branch-specific local sensitivity (Paul15) | `analysis/figure3_hematopoiesis/figure3_hematopoiesis.ipynb` | — |
| **Figure 4** | Concordance vs SpliceJAC and Dynamo (Setty 2019 CD34+); transition-gene overlap under top-HVG null | Four sequential notebooks: `analysis/figure4_concordance/fig4_01_scjdo_operators.ipynb` → `fig4_02_splicejac_operators.ipynb` → `fig4_03_dynamo_operators.ipynb` → `fig4_04_concordance_metrics.ipynb` | — |
| **Figure 5** | Reprogramming SCP295 (IPS, MET, Stromal, IPS_late) with TV=0 default | `analysis/figure5_reprogramming/figure5_reprogramming.ipynb` | — |
| — mild-TV variant (Supp Fig S4 anchor) | Same pipeline with `TV_LAMBDA = 0.1` opt-in | `analysis/figure5_reprogramming/figure5_reprogramming_tv01.ipynb` | — |
| **Figure 6** | K562 CRISPRi Perturb-seq Schrödinger bridges (NT → PVT1, MALAT1, PSMA3-AS1) | `analysis/figure6_k562_perturbseq/figure6_k562.ipynb` | — |

---

## Main-text tables

| Deliverable | Description | Primary code |
|---|---|---|
| **Table 1** | Decision guide (narrative — no computed panels) | *n/a* |

---

## Supplementary figures / notes / tables

| Deliverable | Description | Primary code |
|---|---|---|
| **Supp. Fig. S2** | Cross-validation curves for archetype-count $K$ selection (Fig. 3 and Fig. 5 datasets) | Panels are produced inside `figure3_hematopoiesis.ipynb` and `figure5_reprogramming.ipynb` (search for `cross-val` / `cv_curve` cells) |
| **Supp. Fig. S3** | Bandwidth sweep (DC and MET peak stability across `h`) | Data: `analysis/supp_bandwidth_sweep/bandwidth_sweep_dc_met.py`.  Figure assembly: `analysis/supp_bandwidth_sweep/make_supp_fig_bandwidth.py` |
| **Supp. Fig. S4** | TV=0 vs TV=0.1 archetype-peak preservation on the MET branch | `analysis/figure5_reprogramming/make_supp_fig_tv_robustness.py` (consumes both TV=0 and TV=0.1 runs of `figure5_reprogramming*.ipynb`) |
| **Supp. Note (multiome, cell-cycle)** | 10x Genomics E18 mouse brain multiome, pre- and post-cell-cycle correction; ExcNeuron instability lists | Pre-correction pipeline: `analysis/supp_multiome_note/run_uncorrected_multiome_pipeline.py`.  Notebooks: `analysis/supp_multiome_note/supp_multiome_drift_cellcycle_regressed.ipynb` and `analysis/supp_multiome_note/supp_multiome_FA_integration.ipynb` |
| **Supp. Table 1** | Embedding-benchmark composite (8 methods × 2 specificity rules × timing on/off) | `analysis/figure2_synthetic_and_benchmarks/fa_composite_reanalysis.py`.  *Post-v52 rerun*: `scratchpad_run/fa_composite_reanalysis.py` (identical output; kept in scratchpad for the corrected re-run of Section 3.2) |
| **Supp. Table 2** | Accession numbers (narrative — no computed panels) | *n/a* |

---

## Load-bearing sub-analyses referenced in-text (not their own figure)

| Result | Where cited | Primary code |
|---|---|---|
| Bifurcation-saddle limit (Results §3) — **corrected reruns** are the definitive numbers | Results §3, Discussion, Failure Modes | **Corrected (definitive):** `scratchpad_run/rescore_baseline_corrected.py`, `scratchpad_run/saddle_dense_sampling_corrected.py`, `scratchpad_run/saddle_readouts_midinterval.py`.  **Originals (kept for auditability only; their internal `ALPHA_CRIT = 1.5` does not correspond to the analytic pitchfork):** `analysis/supp_saddle_localization/saddle_readouts.py`, `analysis/supp_saddle_localization/saddle_dense_sampling.py`, `analysis/supp_saddle_localization/saddle_localization_v3_label_blind.ipynb` |
| Palantir–DPT rank-correlation on the erythroid branch ($\rho=0.08$) | Results §Hematopoiesis | `analysis/figure3_hematopoiesis/palantir_dpt_agreement.py` (plus `palantir_dpt_agreement_extra_variants.py` for the variant sweep) |
| Top-HVG-matched null for transition-gene overlap in the concordance analysis (49× fold-over-null) | Results §Concordance | `analysis/figure4_concordance/top_hvg_null.py` (called from `fig4_04_concordance_metrics.ipynb`) |
| TV-robustness controls on synthetic sequential systems (Systems 1a/1b) at `tv_lambda ∈ {0, 0.1, 1, 10}` | Results §Reprogramming (Fig. 5 legend) | `analysis/figure2_synthetic_and_benchmarks/tv_robustness_test.py`, `analysis/figure2_synthetic_and_benchmarks/concurrent_fix_tests.py`, `analysis/figure2_synthetic_and_benchmarks/init_at_gt_probe.py` |
| K562 **matched-perturbation program-content null** on G2M convergence (MALAT1/PSMA3-AS1 vs 19 size-matched non-focal lncRNAs; 171 null pairs) | Results §K562 | `scratchpad_run/k562_perm_null_g2m.py` (result JSON: `k562_perm_null_g2m.json`) |

---

## Full reproduction — order of operations

The main figures are independent of each other except for data-preparation caches. A reasonable run order:

1. `analysis/figure2_synthetic_and_benchmarks/fig2a…fig2f` (synthetic — no external data; ~2–4 h total).
2. `analysis/figure3_hematopoiesis/figure3_hematopoiesis.ipynb` (Paul15 via Scanpy; ~30 min).
3. `analysis/figure4_concordance/fig4_01…fig4_04` in order (Setty 2019 CD34+, needs GSE117498; ~1 h).
4. `analysis/figure5_reprogramming/figure5_reprogramming.ipynb` (SCP295; ~1–2 h). Then rerun with `TV_LAMBDA = 0.1` or open the pre-built `figure5_reprogramming_tv01.ipynb`.
5. `analysis/figure6_k562_perturbseq/figure6_k562.ipynb` (Replogle 2022 K562 lncRNA sublibrary; ~30 min).
6. Supplementary: `analysis/supp_multiome_note/*`, `analysis/supp_bandwidth_sweep/*`.
7. **After the main figures**, the corrected saddle reruns and the K562 null in `scratchpad_run/`. Each is ~10–15 min:
   - `scratchpad_run/rescore_baseline_corrected.py` (rescores saved curves — <1 s)
   - `scratchpad_run/saddle_dense_sampling_corrected.py` (2 seeds ~30 min)
   - `scratchpad_run/saddle_readouts_midinterval.py` (3 seeds ~30 min)
   - `scratchpad_run/k562_perm_null_g2m.py` (19 bridges ~7 min)

`analysis/run_all.sh` is a convenience wrapper for the analysis/ subtree; it does **not** invoke the scratchpad_run corrections and must be complemented by them for the v52 numbers.

---

## Data dependencies

| Dataset | Access | Used by |
|---|---|---|
| Paul15 hematopoiesis | `scanpy.datasets.paul15()` | Fig. 2 embedding benchmark, Fig. 3, palantir_dpt_agreement |
| Setty 2019 CD34+ bone marrow | GEO **GSE117498** | Fig. 4 concordance |
| SCP295 iPSC reprogramming (Schiebinger 2019) | Single Cell Portal SCP295 | Fig. 5 |
| Replogle 2022 K562 CRISPRi Perturb-seq (lncRNA sublibrary) | See paper's Data Availability | Fig. 6, K562 null |
| 10x Genomics E18 mouse brain 5k multiome | 10xgenomics.com | Supp. Note (multiome cell-cycle) |
| Synthetic SDE systems | Generated in-script | Fig. 2, all saddle-localization scripts |

`analysis/cache/` and `analysis/results/` are locally generated and are gitignored — they will be recreated by the notebooks on first run.
