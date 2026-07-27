# Figure paths — where every manuscript figure lives

Quick-reference mapping from **manuscript figure number → file paths → source code**.
All paths are relative to the repo root (`/Users/terooatt/Downloads/scJDO/`).
Every figure exists as both `.pdf` (vector) and `.png` (300 DPI raster).

## Main figures

| Paper | PDF | PNG | Source code | Regenerate command (from repo root) |
|---|---|---|---|---|
| **Fig 1** | — schematic, hand-drawn, no code — | | | |
| **Fig 2** | `final_manuscript/figures/main/Figure2.pdf` | `final_manuscript/figures/main/Figure2.png` | `final_manuscript/code/figure2_synthetic_and_benchmarks/fig2b_peak_timing_null_analysis.ipynb` | `jupyter nbconvert --execute --to notebook <path>` |
| **Fig 3** | `final_manuscript/figures/main/Figure3.pdf` | `final_manuscript/figures/main/Figure3.png` | `final_manuscript/code/figure3_hematopoiesis/figure3_hematopoiesis.ipynb` | `jupyter nbconvert --execute --to notebook <path>` |
| **Fig 4** | `final_manuscript/figures/main/Figure4.pdf` | `final_manuscript/figures/main/Figure4.png` | `final_manuscript/code/figure4_concordance/fig4_04_concordance_metrics.ipynb` | `jupyter nbconvert --execute --to notebook <path>` |
| **Fig 5** (TV=0, main) | `final_manuscript/figures/main/Figure5_TV0.pdf` | `final_manuscript/figures/main/Figure5_TV0.png` | `final_manuscript/code/figure5_reprogramming/figure5_reprogramming.ipynb` | `jupyter nbconvert --execute --to notebook <path>` |
| **Fig 5** (TV=0.1, companion) | `final_manuscript/figures/main/Figure5_TV01.pdf` | `final_manuscript/figures/main/Figure5_TV01.png` | `final_manuscript/code/figure5_reprogramming/figure5_reprogramming_tv01.ipynb` | `jupyter nbconvert --execute --to notebook <path>` |
| **Fig 6** (main) | `final_manuscript/figures/main/Figure6.pdf` | `final_manuscript/figures/main/Figure6.png` | `final_manuscript/code/figure6_k562_perturbseq/figure6_k562.ipynb` | `jupyter nbconvert --execute --to notebook <path>` |
| **Fig 6** (validation panel) | `final_manuscript/figures/main/Figure6_validation.pdf` | `final_manuscript/figures/main/Figure6_validation.png` | (same notebook as Fig 6 above) | |

## Supplementary figures

| Description | PDF | PNG | Source code |
|---|---|---|---|
| Multiome joint FA integration | `final_manuscript/figures/supplementary/FigureS_multiome_FA_integration.pdf` | `.png` | `final_manuscript/code/supp_multiome_note/supp_multiome_FA_integration.ipynb` |
| Multiome integration QC panels | `final_manuscript/figures/supplementary/FigureS_multiome_integration_QC.pdf` | `.png` | (same notebook as above) |
| Multiome marker ordering pre-vs-post cell-cycle regression | `final_manuscript/figures/supplementary/FigureS_multiome_marker_ordering_pre_vs_post_cc.pdf` | `.png` | `final_manuscript/code/supp_multiome_note/supp_multiome_drift_cellcycle_regressed.ipynb` |
| Multiome ExcNeuron instability, cell-cycle-regressed | `final_manuscript/figures/supplementary/FigureS_multiome_instab_ExcNeuron_cc_regressed.pdf` | `.png` | (same notebook as above) |
| Multiome ExcNeuron instability, uncorrected | `final_manuscript/figures/supplementary/FigureS_multiome_instab_ExcNeuron_uncorrected.pdf` | `.png` | `final_manuscript/code/supp_multiome_note/run_uncorrected_multiome_pipeline.py` |
| Bandwidth sweep — DC / MET | `final_manuscript/figures/supplementary/FigureS_bandwidth_sweep_DC_MET.pdf` | `.png` | `final_manuscript/code/supp_bandwidth_sweep/bandwidth_sweep_dc_met.py` + `make_supp_fig_bandwidth.py` |
| Saddle-localization label-blind test | `final_manuscript/figures/supplementary/FigureS_saddle_localization_label_blind.pdf` | `.png` | `final_manuscript/code/supp_saddle_localization/saddle_localization_v3_label_blind.ipynb` |

## Individual panels + diagnostics

Located in `final_manuscript/figures/panels_and_diagnostics/` — 80 files (40 × PDF + PNG).
Use these if you need to redo a single sub-panel of a composite main figure, or want to inspect intermediate diagnostic outputs.

Naming: `fig<N>_<subject>.{pdf,png}` — e.g., `fig3_instab_Ery.pdf`, `fig6_MALAT1_reg_summary_forward.pdf`.
See `final_manuscript/figures/README.md` for the full per-file description.

## How to (re)generate any figure

### Option 1 — regenerate a single figure

The PDF file is produced by the notebook / script listed in the tables above. From the repo root:

```bash
# Notebook (e.g., Fig 3)
cd final_manuscript/code
jupyter nbconvert --to notebook --execute --inplace figure3_hematopoiesis/figure3_hematopoiesis.ipynb

# Python script (e.g., Fig 5 bandwidth-sweep supp)
cd final_manuscript/code/supp_bandwidth_sweep
python bandwidth_sweep_dc_met.py
python make_supp_fig_bandwidth.py
```

The output PDF lands in `analysis/results/figure<N>/` (or a sibling location — see the notebook's `OUTDIR = …` line). Copy it into `final_manuscript/figures/main/` (or `supplementary/`) under the paper-figure name from the table above.

### Option 2 — regenerate the full chain

```bash
cd final_manuscript/code
./run_all.sh              # ~8 h; runs all 27 items sequentially
./run_all.sh fast         # ~2 h; skips the three heaviest notebooks
./run_all.sh fig4         # just Fig 4 items (~8 min)
```

Per-item logs land in `final_manuscript/code/logs/`.

### Regenerate PNG from PDF (300 DPI)

If you've updated a PDF and need a matching PNG:

```bash
pdftoppm -r 300 -png -singlefile final_manuscript/figures/main/Figure3.pdf \
                                 final_manuscript/figures/main/Figure3
```

The `-singlefile` flag ensures the output is `Figure3.png` (not `Figure3-1.png`).

## Data required to regenerate

Fig 3 and Fig 4 use the shipped marrow data (`examples/marrow_sample_scseq_counts.h5ad`) — no extra setup needed.

Fig 5 (SCP295), Fig 6 (K562 Perturb-seq), and the Supplementary Note multiome figures need external datasets — see `final_manuscript/DATA.md` for download links and either the exact paths to place them or the environment variables (`SCJDO_SCP295_H5AD`, `SCJDO_K562_DATA_ROOT`, `SCJDO_MULTIOME_MTX_DIR`) that redirect the notebooks to your copy.

## Verification — do these figures match the manuscript numbers?

Yes. Every quantitative claim in `manuscript/scJDO_manuscript_v56.md` is pinned to the corresponding notebook — see the **"Manuscript-to-code number pins"** table in `final_manuscript/README.md`. All main figures were regenerated in the 2026-07-24 full re-run; Fig 3 instab panels were regenerated on 2026-07-26 after the auto-relax threshold fix in `scjdo/pl/_drift.py`.
