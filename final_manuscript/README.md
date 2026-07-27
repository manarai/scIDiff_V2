# scJDO — Reviewer replication package

This directory contains the final manuscript (v56) and a self-contained
copy of the analysis code that generates every figure. The chain has been
end-to-end re-executed against manuscript v56 (2026-07-24 → 2026-07-26)
and every cited number matches the reproducible output.

## Contents

    final_manuscript/
      README.md                 (this file)
      DATA.md                   (data availability + setup)
      environment.yml           (conda environment spec)
      manuscript/
        scJDO_manuscript_v56.md
      figures/                  (all main + supplementary figures, PDF + 300 DPI PNG)
        README.md               (per-figure source-notebook mapping)
        main/                   Figure2 … Figure6 (+ Fig5 TV=0.1 companion, Fig6 validation)
        supplementary/          Supplementary Note + bandwidth-sweep + saddle-localization
        panels_and_diagnostics/ Individual panels + diagnostic subfigures
      code/
        run_all.sh              (sequential runner for the full chain)
        examples -> ../../examples/   (symlink; must resolve — see DATA.md)
        figure2_synthetic_and_benchmarks/
        figure3_hematopoiesis/
        figure4_concordance/
        figure5_reprogramming/
        figure6_k562_perturbseq/
        supp_bandwidth_sweep/
        supp_multiome_note/
        supp_saddle_localization/

## Figure → notebook chain (paper-figure numbering)

| Paper figure  | Producer(s)                                                                              |
| ------------- | ---------------------------------------------------------------------------------------- |
| Fig 1         | schematic — no notebook                                                                  |
| Fig 2         | `figure2_synthetic_and_benchmarks/fig2a…fig2f_*.ipynb` + supporting `.py` scripts       |
| Fig 3         | `figure3_hematopoiesis/figure3_hematopoiesis.ipynb` + `palantir_dpt_agreement*.py`      |
| Fig 4         | `figure4_concordance/fig4_01→fig4_04_*.ipynb` + `top_hvg_null.py`                       |
| Fig 5         | `figure5_reprogramming/figure5_reprogramming.ipynb` (TV=0) + `_tv01.ipynb` (TV=0.1)     |
| Fig 6         | `figure6_k562_perturbseq/figure6_k562.ipynb`                                             |
| Supp. Note    | `supp_multiome_note/supp_multiome_FA_integration.ipynb` + `_drift_cellcycle_regressed.ipynb` |
| Supp. figs    | `supp_bandwidth_sweep/`, `supp_saddle_localization/`                                     |

## Quick start (one command)

Prerequisites (see `DATA.md` and `environment.yml`):

1. Install `scjdo` (v0.3.0) and dependencies — `pip install -e ..` from the
   parent scJDO repo, or follow `environment.yml`.
2. Ensure the marrow data file `../examples/marrow_sample_scseq_counts.h5ad`
   exists (or point the `examples` symlink at your copy).
3. Set `SCP295_H5AD`, `K562_TAR_DIR`, and `MULTIOME_MTX_DIR` for the
   dataset-specific paths (see `DATA.md`).

Then:

    cd code/
    ./run_all.sh              # full run (~8 h wall-clock on M-series Mac)
    ./run_all.sh fast         # skip the 3 heaviest notebooks (~2 h)
    ./run_all.sh fig4         # run only Fig 4 items

Status lands in `code/run_all.status`; per-item logs in `code/logs/`.
Executed notebooks are rewritten in place with fresh outputs.

## Manuscript-to-code number pins

Every quantitative claim in v56 can be regenerated from this tree:

| Manuscript number                                    | Regenerate by                                        |
| ---------------------------------------------------- | ---------------------------------------------------- |
| Fig 2a S1/S3/S4 (0.0005 / 0.228 / 0.065)             | `figure2_synthetic_and_benchmarks/fig2a_synthetic_benchmark.ipynb` |
| Fig 2 concurrent-prev 0.025                          | `figure2_synthetic_and_benchmarks/concurrent_fix_tests.py` (K=5 row) |
| Fig 2 GT-init residual 0.1101 vs 0.1105              | `figure2_synthetic_and_benchmarks/init_at_gt_probe.py` |
| Fig 2 TV=0.1 preserves seq. τ=1.00 (8/8)             | `figure2_synthetic_and_benchmarks/tv_robustness_test.py` |
| Saddle Re(λ) 0.83, P⊥JP⊥ 0.85, FTLE 0.78, div ±0.05 | `supp_saddle_localization/saddle_readouts.py`        |
| Saddle dense-enrichment 0.93 / 0.86 / 0.86           | `supp_saddle_localization/saddle_dense_sampling.py`  |
| Fig 3 branch counts + R² + GATA1                     | `figure3_hematopoiesis/figure3_hematopoiesis.ipynb` |
| Fig 3 Palantir/DPT Ery ρ ∈ [0.839, 0.988]            | `figure3_hematopoiesis/palantir_dpt_agreement*.py`  |
| Fig 4 Jaccard 0.0749 / 0.0051; fold 32× / 2.2×       | `figure4_concordance/fig4_04_concordance_metrics.ipynb` |
| Fig 4 top-HVG-null fold 10× / 0.7×                   | `figure4_concordance/top_hvg_null.py` (reads fig4_04's JSON) |
| Fig 4 ρ = -0.818 (Precursors argmax)                 | `figure4_concordance/fig4_04_concordance_metrics.ipynb` |
| Fig 5 MET τ_peak 0.77; TV=0/0.1 archetype peaks     | `figure5_reprogramming/figure5_reprogramming{,_tv01}.ipynb` |
| Fig 5 bandwidth sweep h=0.01→0.15                    | `supp_bandwidth_sweep/bandwidth_sweep_dc_met.py`    |
| Fig 6 MALAT1 FDR 2e-17 / 12 mitotic genes            | `figure6_k562_perturbseq/figure6_k562.ipynb`         |
| Fig 6 PVT1 permutation p=0.143                       | `figure6_k562_perturbseq/figure6_k562.ipynb`         |

## Reference outputs

`figure4_concordance/concordance_metrics.json` is checked in — it is the
frozen Fig 4 output that the manuscript's Fig 4 numbers cite verbatim.
`top_hvg_null.py` reads from this file, so if a reviewer wants the two-null
fold decomposition **without** rerunning the full pipeline, they can run
`top_hvg_null.py` on the checked-in JSON directly and reproduce the
`Dynamo 32× / 10×; SpliceJAC 2.2× / 0.7×` line.

## Wall-clock estimates (per-item, single machine)

The heavy tail:

- `fig2f_fa_pca_scvi_benchmark.ipynb`   ~95 min
- `figure5_reprogramming.ipynb`         ~50 min
- `figure5_reprogramming_tv01.ipynb`    ~50 min
- `fig2e_fa_benchmark.ipynb`            ~80 min
- `fig2d_pseudotime_embedding_benchmark.ipynb`   ~65 min

Everything else is < 20 min per item; total ~8 h for `run_all.sh` on an
M-series Mac (measured 2026-07-24).

## Contact

Corresponding author: see manuscript.
