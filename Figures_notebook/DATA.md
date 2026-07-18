# Data availability for `Figures_notebook/`

This document lists every dataset the figure notebooks depend on, where each
was obtained, the expected on-disk path used by the current notebooks, and any
preprocessing required to reach the input files each notebook actually reads.

Reviewers reproducing the figures should either recreate the paths below, or
edit the `DATA_PATH` / `DATA_ROOT` constants at the top of each notebook to
point at a local copy.

---

## Figure 3 — mouse bone-marrow hematopoiesis

**Dataset.** The 4,142-cell mouse bone-marrow sample distributed with the
Palantir tutorial (from the Setty *et al.* 2019 CD34+ profiling; also included
as scanpy example data). **This is not** `sc.datasets.paul15()` — the Paul et
al. 2015 dataset would have a different cell count and a different set of
progenitor labels.

**Source.** `https://github.com/dpeerlab/Palantir` — `data/marrow_sample_scseq_counts.h5ad`
(shipped in the Palantir repository under `data/`).

**Expected local path.**
```
/Users/terooatt/Downloads/scJDO/examples/marrow_sample_scseq_counts.h5ad
```
Loaded by `Figure3_FA.ipynb` cell 4:
```python
DATA_PATH = '../examples/marrow_sample_scseq_counts.h5ad'
```

**Terminal-state anchors used.**
- Root cell (progenitor): `Run5_164698952452459`
- Terminal DC:            `Run5_131097901611291`
- Terminal Mono:          `Run5_134936662236454`
- Terminal Ery:           `Run4_200562869397916`

**Post-Palantir branch cell counts (Ery / DC / Mono).** 1,151 / 1,903 / 2,008
(reproduces exactly under `q=0.01, eps=0.01` in `palantir.presults.select_branch_cells`).

---

## Figure 4 (concordance analysis) — Setty CD34+ bone marrow

**Dataset.** Setty *et al.* 2019 CD34+ human bone marrow with spliced/unspliced
layers, distributed via `scvelo.datasets.bonemarrow()`. 5,780 cells after
`scv.pp.filter_and_normalize(min_shared_counts=20)`, 6,482 measured genes.

**Source.** Downloaded on first call by:
```python
import scvelo as scv
ad = scv.datasets.bonemarrow()          # scVelo caches under ~/.cache/scvelo/
```

**No manual download required** — the concordance notebooks
(`Manuscript/concordance_01_scjdo.ipynb`, `_02_splicejac.ipynb`,
`_03_dynamo.ipynb`) all resolve it via scvelo.

---

## Figure 5 — SCP295 iPSC reprogramming (Schiebinger 2019)

**Dataset.** SCP295 Single-Cell Portal deposit accompanying
Schiebinger *et al.* 2019 (Cell). Dense MEF → iPSC time course
(day 0 to day 18, ~0.5-day interval). The notebook uses a **24,999-cell
stratified subsample** — this is the shipped input, not the full raw upload.

**Source.** Broad Single-Cell Portal:
`https://singlecell.broadinstitute.org/single_cell/study/SCP295`
(requires login; raw expression + `cell_sets.gmt` + `2i.h5ad` / `serum.h5ad`).

**Preprocessing to the notebook's input file.** Run
`Figures_notebook/build_scp295_h5ad.py`, which pools the day-wise files,
filters to the four annotated branches (`IPS`, `MET`, `Stromal`, `IPS_late`),
computes the day-normalized `day_norm ∈ [0, 1]`, and builds the FLE embedding.
The output is a single `scp295.h5ad` (24,999 cells × 19,089 genes; obs includes
`day`, `day_norm`, `cell_set`; obsm includes `X_fle`).

**Expected local path** (as currently hard-coded at the top of Figure4.ipynb / Figure4_tv01.ipynb):
```
/Users/terooatt/Documents/Project_scQDiff/02_scQDiff/scIDIFF_anndata/data/SCP295/scp295.h5ad
```

**Cell-count reconciliation note.** The manuscript v38 currently states
n = 27,538 for this dataset; the shipped `scp295.h5ad` is n = 24,999 in the
notebook prose, in the panel-a title, and in `results/figure4/branch_summary.csv`.
The authoritative number is **24,999** (from `build_scp295_h5ad.py`).

---

## Figure 6 — K562 CRISPRi Perturb-seq (10x Genomics)

**Dataset.** 10x Genomics 84K K562 5'v3 GEM-X CRISPR aggregate
(KRAB-dCas9 K562 cells transduced with a CRISPRi lncRNA + protein-coding
sgRNA library, ~6,903 sgRNAs). Non-targeting (NT) control + three
lncRNA-knockdown populations used (PVT1, MALAT1, PSMA3-AS1).

**Source.** 10x Genomics dataset page (CC BY 4.0):
`https://www.10xgenomics.com/datasets/84k_K562_5pv3_GEMX_CRISPR_aggregate`

Two tarballs required:
1. `filtered_feature_bc_matrix.tar.gz` (~1.2 GB) — GEX + sgRNA feature matrix
2. `crispr_analysis.tar.gz` — Cell Ranger CRISPR call table

**Expected local root.**
```
/Users/terooatt/Documents/Project_scQDiff/02_scQDiff/scIDIFF_anndata/data/pertubeseq_10x/
    ├── <filtered_feature_bc_matrix.tar.gz>
    └── <crispr_analysis.tar.gz>
```
`Figure5.ipynb` cells 3–4 idempotently extract both tarballs to an `extracted/`
subdirectory and cache the loaded AnnData under `Figures_notebook/cache/`.

---

## Supplementary Note — 10x E18 mouse-brain multiome

**Dataset.** 10x Genomics E18 mouse brain fresh 5k Multiome. Joint RNA + ATAC.

**Source.** 10x Genomics dataset page:
`https://www.10xgenomics.com/datasets/fresh-embryonic-e-18-mouse-brain-5-k-1-standard-2-0-0`

Requires the `filtered_feature_bc_matrix/` directory (Cell Ranger ARC output).

**Expected local root.**
```
/Users/terooatt/Documents/Project_scQDiff/02_scQDiff/scIDIFF_anndata/data/multiomic/filtered_feature_bc_matrix/
```
Loaded by `Figure6_02_multiome_FA_integration.ipynb` cell 4 via
`sc.read_10x_mtx(..., gex_only=False)`.

**Downstream pipeline.**
1. `Figure6_02_multiome_FA_integration.ipynb` — reads the raw 10x matrix,
   splits GEX / ATAC, computes per-modality FA (RNA and ATAC), joins them,
   and writes:
   ```
   Figures_notebook/results/figure6_multiome_fa/adata_multiome_fa.h5ad
   Figures_notebook/results/figure6_multiome_fa/adata_multiome_fa_atac.h5ad
   ```
2. `Figure6_01_multiome_drift_cellcycle_regressed_FA.ipynb` — reads those two
   h5ads and fits the CC-regressed drift/instability pipeline.

**Cell-cycle-correction reproducibility note.** The joint FA input shipped in
`adata_multiome_fa.h5ad` already has per-modality S/G2M regression applied
upstream (all three FA latents have Pearson |r| ≈ 0 with S/G2M scores at the
leading components). A *truly-uncorrected* pipeline for the multiome
Supplementary Note requires re-running `Figure6_02` with S/G2M regression
disabled in the per-modality FA step; this reproduction is documented in
`analysis/results/supp_multiome/precorrection_reproducibility_check.json`.

---

## Where the figures go once regenerated

Every `results/` subdirectory under `Figures_notebook/` is committed with the
current-run PDFs, PNGs, CSVs, and JSONs. If you regenerate, please version any
new outputs into a per-figure subdirectory rather than overwriting shared
files.

---

## Dataset ↔ figure quick reference

| Figure | Dataset | External path (root)                                                       | Loader |
|-------:|---------|----------------------------------------------------------------------------|--------|
| 2      | synthetic (no external data)             | —                                                          | in-notebook simulation |
| 3      | mouse marrow sample (Palantir tutorial)  | `scJDO/examples/`                                          | `sc.read` from repo    |
| 4      | Setty CD34+ (scvelo bonemarrow)          | `~/.cache/scvelo/` (auto-download)                         | `scvelo.datasets.bonemarrow` |
| 5      | SCP295 iPSC reprogramming                | `~/Documents/…/scIDIFF_anndata/data/SCP295/`               | `sc.read` cached h5ad  |
| 6      | K562 84K CRISPRi Perturb-seq             | `~/Documents/…/scIDIFF_anndata/data/pertubeseq_10x/`       | `sc.read_10x_mtx` + tar extract |
| S (multiome) | 10x E18 mouse brain multiome       | `~/Documents/…/scIDIFF_anndata/data/multiomic/`            | `sc.read_10x_mtx` |
