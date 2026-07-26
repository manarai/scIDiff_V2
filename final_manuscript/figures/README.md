# Figures — final_manuscript/figures/

All figures for scJDO manuscript v54, at high resolution.

Every file exists as both `.pdf` (vector) and `.png` (300 DPI raster).
The PNGs were generated from the PDFs with `pdftoppm -r 300 -png -singlefile`;
regenerate any of them with the same command if the PDF changes.

## Directory layout

    figures/
      README.md                       (this file)
      main/                           7 main paper figures (14 files: PDF + PNG)
      supplementary/                  Supplementary Note & supplementary figures
      panels_and_diagnostics/         Individual panels + diagnostic subfigures

## Main figures

Fig 1 is a hand-drawn conceptual schematic; it has no code-generated version.
All other main figures are code-generated and pinned below.

| File | Manuscript | Producer notebook | Notes |
|---|---|---|---|
| `main/Figure2.pdf` / `.png` | Fig 2 (composite) | `code/figure2_synthetic_and_benchmarks/fig2b_peak_timing_null_analysis.ipynb` | Manuscript-formatted composite of Fig 2 panels |
| `main/Figure3.pdf` / `.png` | Fig 3 | `code/figure3_hematopoiesis/figure3_hematopoiesis.ipynb` | Hematopoiesis (Palantir/Setty marrow) |
| `main/Figure4.pdf` / `.png` | Fig 4 | `code/figure4_concordance/fig4_04_concordance_metrics.ipynb` | Dynamo / SpliceJAC concordance |
| `main/Figure5_TV0.pdf` / `.png` | Fig 5 (TV=0, main) | `code/figure5_reprogramming/figure5_reprogramming.ipynb` | SCP295 reprogramming, unregularised semi-NMF |
| `main/Figure5_TV01.pdf` / `.png` | Fig 5 companion (TV=0.1) | `code/figure5_reprogramming/figure5_reprogramming_tv01.ipynb` | Same pipeline with `tv_lambda = 0.1` |
| `main/Figure6.pdf` / `.png` | Fig 6 (main) | `code/figure6_k562_perturbseq/figure6_k562.ipynb` | K562 CRISPRi Perturb-seq bridges |
| `main/Figure6_validation.pdf` / `.png` | Fig 6 validation | `code/figure6_k562_perturbseq/figure6_k562.ipynb` | Permutation null + G2M enrichment panels |

## Supplementary figures

| File | Description | Producer |
|---|---|---|
| `supplementary/FigureS_multiome_FA_integration.pdf` / `.png` | Joint RNA+ATAC FA integration | `code/supp_multiome_note/supp_multiome_FA_integration.ipynb` |
| `supplementary/FigureS_multiome_integration_QC.pdf` / `.png` | Multiome integration QC panel | (same as above) |
| `supplementary/FigureS_multiome_marker_ordering_pre_vs_post_cc.pdf` / `.png` | Cell-cycle correction: marker-ordering diff | `code/supp_multiome_note/supp_multiome_drift_cellcycle_regressed.ipynb` |
| `supplementary/FigureS_multiome_instab_ExcNeuron_cc_regressed.pdf` / `.png` | ExcNeuron instability after CC regression | (same as above) |
| `supplementary/FigureS_multiome_instab_ExcNeuron_uncorrected.pdf` / `.png` | ExcNeuron instability without CC regression | `code/supp_multiome_note/run_uncorrected_multiome_pipeline.py` |
| `supplementary/FigureS_bandwidth_sweep_DC_MET.pdf` / `.png` | Bandwidth sensitivity (DC / MET), pinned in Fig 5 discussion | `code/supp_bandwidth_sweep/bandwidth_sweep_dc_met.py` + `make_supp_fig_bandwidth.py` |
| `supplementary/FigureS_saddle_localization_label_blind.pdf` / `.png` | Label-blind flow-topology saddle localization test | `code/supp_saddle_localization/saddle_localization_v3_label_blind.ipynb` |

## Panels and diagnostics (`panels_and_diagnostics/`)

Individual panels and intermediate diagnostic subfigures. Useful if a
reviewer wants to inspect a single panel from a composite figure, or wants
to see the diagnostic outputs of individual synthetic-benchmark stages.
Not intended for inclusion in the manuscript body — the main-figure files
above already composite these appropriately.

**Fig 2 diagnostics** (14 files, `fig2_*`):

- `fig2_fig01_true_lambda`, `fig2_fig02_sweep_A_lambda_corr`, … through
  `fig2_fig14a_bifurcation_v3_ground_truth` — the individual per-stage
  outputs of `fig2a_synthetic_benchmark.ipynb`.
- `fig2_synthetic_benchmark_composite` — the raw 4-panel composite before
  layout for the paper.
- `fig2_pipeline_recovery_probe`, `fig2_synthetic_peak_timing_system_overview`,
  `fig2_synthetic_peak_timing_target_pvalues` — panels supporting the
  concurrent-activation and peak-timing-null discussion in the Results.
- `fig2c_embedding_benchmark_8method`, `fig2e_fa_benchmark`,
  `fig2f_fa_pca_scvi_benchmark` — the individual embedding-comparison
  figures from fig2c/e/f.

**Fig 3 panels** (5 files, `fig3_*`):

- `fig3_instab_{Ery,DC,Mono}` — per-branch instability-gene panels.
- `fig3_{Ery,DC}_reg_summary` — per-branch regulator summary panels.

**Fig 5 panels** (8 files, `fig5*`):

- `fig5_instab_{IPS,IPS_late,MET,Stromal}` — per-branch instability panels
  (TV=0 pipeline).
- `fig5_tv01_instab_{IPS,IPS_late,MET,Stromal}` — same, TV=0.1 pipeline.

**Fig 6 panels** (5 files, `fig6_*`):

- `fig6_umap_perturbations` — the joint UMAP of NT + three perturbation
  populations.
- `fig6_{PVT1,MALAT1,PSMA3-AS1}_reg_summary_forward` — per-perturbation
  regulator summaries for the forward bridge (NT → target).

## Regenerating a figure

To regenerate any single figure:

1. Ensure the corresponding notebook / script has run (see the pin-table in
   `code/…` or `README.md`).
2. Find the source PDF in the notebook's output directory (typically
   `analysis/results/figureN/` when running from a full scJDO clone).
3. Copy the PDF into `figures/main/` (or `supplementary/`, or
   `panels_and_diagnostics/`) under the appropriate name.
4. Generate the 300 DPI PNG:

        pdftoppm -r 300 -png -singlefile figures/main/FigureN.pdf figures/main/FigureN

The full chain regenerates every figure in one shot:

    cd ../code
    ./run_all.sh            # ~8 h on M-series Mac
    # then re-run the copy step for any updated PDFs
