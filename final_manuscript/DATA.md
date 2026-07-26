# Data availability and reviewer setup

Every dataset used by the code chain is public. Two datasets download
automatically at runtime; the other three (SCP295, K562 Perturb-seq, and
the multiome) must be placed at a location the notebooks can read.

## Datasets used

| Dataset                    | Source / accession                                                                                          | Used by                                                             | Auto-download? |
| -------------------------- | ------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------- | -------------- |
| Palantir marrow (4,142 cells) | included as `examples/marrow_sample_scseq_counts.h5ad` in the scJDO repo (originally Setty *et al.* 2019)   | Fig 2c/2d/2e/2f, Fig 3                                              | shipped with repo |
| Setty 2019 CD34+ bone marrow (5,780 cells)  | `scvelo.datasets.bonemarrow()` — auto-download via scvelo                                                    | Fig 4 concordance chain                                             | yes (first run) |
| SCP295 iPSC reprogramming (24,999 cells, days 0–18) | Schiebinger *et al.* 2019, SCP295 on Single Cell Portal — H5AD conversion required (see below) | Fig 5 (both TV notebooks), `supp_bandwidth_sweep/bandwidth_sweep_dc_met.py` | no |
| K562 CRISPRi Perturb-seq (83,943 cells)  | 10x Genomics `84K_K562_5pv3_GEMX_CRISPR_aggregate` — three `.tar.gz` files (CC BY 4.0)                       | Fig 6                                                                | no |
| Chen 2019 multiome (10x Multiome, 10k cells) | 10x Genomics filtered feature-barcode matrix                                                                 | Supplementary Note (both notebooks)                                  | no |

## Where each dataset must live

The code assumes these paths. **Two options** for each: either place the
file at the exact path below, or edit the path at the top of the notebook
before running.

### 1. Marrow (Setty 2019, Palantir tutorial)

    <repo-root>/examples/marrow_sample_scseq_counts.h5ad

Notebooks reference this via `../examples/marrow_sample_scseq_counts.h5ad`,
which resolves through the `code/examples` symlink to
`<repo-root>/examples/`. **If you moved this directory,** re-point the
symlink:

    cd final_manuscript/code
    rm examples && ln -s /path/to/your/examples examples

### 2. Setty 2019 CD34+ (Fig 4)

Downloaded automatically on first execution of `fig4_04_concordance_metrics.ipynb`
via `scvelo.datasets.bonemarrow()`. Cached under `~/.scvelo/`. No action
required.

### 3. SCP295 (Fig 5)

Fig 5 notebooks read:

    DATA_PATH = '/Users/terooatt/Documents/Project_scQDiff/02_scQDiff/scIDIFF_anndata/data/SCP295/scp295.h5ad'

You **must** either:
- Download SCP295 from the Broad Single Cell Portal, convert to `.h5ad`
  (24,999 cells × 19,089 genes, `obs['day']` and `obs['cell_set']` present),
  and place at the exact path above, **or**
- Edit `DATA_PATH` at the top of `figure5_reprogramming/figure5_reprogramming.ipynb`
  and `figure5_reprogramming_tv01.ipynb` and
  `supp_bandwidth_sweep/bandwidth_sweep_dc_met.py`.

### 4. K562 Perturb-seq (Fig 6)

`figure6_k562.ipynb` reads:

    DATA_ROOT   = '/Users/terooatt/Documents/Project_scQDiff/02_scQDiff/scIDIFF_anndata/data/pertubeseq_10x'
    TAR_MTX     = <DATA_ROOT>/84K_K562_5pv3_GEMX_CRISPR_aggregate_count_filtered_feature_bc_matrix.tar.gz
    TAR_CRISPR  = <DATA_ROOT>/84K_K562_5pv3_GEMX_CRISPR_aggregate_count_crispr_analysis.tar.gz

Download the three `.tar.gz` files from
`https://www.10xgenomics.com/datasets/84k_K562_5pv3_GEMX_CRISPR_aggregate`
into a directory of your choice, then either place at the above path or
edit `DATA_ROOT` at the top of the notebook.

### 5. Multiome (Supplementary Note)

`supp_multiome_note/supp_multiome_FA_integration.ipynb` reads:

    DATA_DIR = '/Users/terooatt/Documents/Project_scQDiff/02_scQDiff/scIDIFF_anndata/data/multiomic/filtered_feature_bc_matrix'

Point this at any 10x Multiome `filtered_feature_bc_matrix` directory
(matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz) that includes both
Gene Expression and Peaks feature types. The Supp. Note pipeline is
tolerant of dataset choice; the shipped analysis used a Chen 2019
10k-cell multiome sample.

## First-time setup checklist

1. Clone the scJDO repo and `pip install -e .` (or follow `environment.yml`).
2. Confirm `examples/marrow_sample_scseq_counts.h5ad` exists in the repo.
3. Verify the symlink:

        ls -la final_manuscript/code/examples
        # should print:  examples -> ../../examples/

4. Download and place SCP295, K562, and multiome data per §3–§5 above
   (or edit the three notebooks' paths).
5. Run `final_manuscript/code/run_all.sh fast` for a first check (~2 h,
   skips the three heaviest notebooks). If everything OKs, run the full
   chain: `run_all.sh` (~8 h).

## What "reproducible" means here

- All 27 items in `run_all.sh` ran end-to-end on 2026-07-24 (macOS 25.5.0,
  Python 3.13.2, scjdo v0.3.0).
- Every quantitative claim in the manuscript that cites a specific number
  is regenerated by one of these items; see the pin-table in `README.md`.
- Where a claim depends on stochastic training (Fig 5 archetype
  decomposition; multiome FA integration), fixed random seeds are used
  and archetype labels are matched by cosine similarity on the operator
  matrices $H$ (see manuscript §Methods and Fig 2 legend).
- The Fig 4 concordance numbers in the manuscript are checked in as
  `code/figure4_concordance/concordance_metrics.json`; reviewers who want
  to verify without a full pipeline run can just `python top_hvg_null.py`
  against that file.
