# Example Notebooks

Five end-to-end tutorials covering all scJDO workflows.
Each notebook uses the high-level API — 4–6 function calls, no boilerplate.

| # | Notebook | Dataset | Analysis | Key outputs |
|---|---|---|---|---|
| 01 | `01_paul15_hybrid_drift_tutorial` | Paul15 hematopoiesis | Drift field + archetypes + instability genes + regulators | `summary_figure`, instability CSV, regulator network |
| 02 | `02_paul15_fourier_tutorial` | Paul15 hematopoiesis | Fourier-domain score network, spectral band decomposition, generated cell validation | Real vs generated power spectrum, band gene lists |
| 03 | `03_schrodinger_bridge_tutorial` | 2D synthetic | Bridge on Gaussian source → target, forward/backward instability, de novo network edges | Trajectory plots, instability asymmetry table |
| 04 | `04_paul15_schrodinger_tutorial` | Paul15 hematopoiesis | Bridge in 50-D PCA space, forward (differentiation) vs backward (de-differentiation) instability genes and regulators | Gene tables (forward/backward), hybrid regulator network |
| 05 | `05_scopatlas_complete_workflow` | Your data + pre-trained model | Stable-operator atlas: operator metrics, embedding, clustering, biological interpretation | Operator regime UMAP, clustering comparison, regime statistics |
| 10 | `10_paul15_koopman_tutorial` | Paul15 hematopoiesis | Koopman spectral backend on the Jacobian tensor: eigenvalues, growth rates, oscillation frequencies, side-by-side comparison with semi-NMF | Complex-plane spectrum, per-mode growth curves, local vs global operator, activation panels |

---

## Running the notebooks

```bash
conda activate scJDO
jupyter lab
```

Start with **01** — it covers the core workflow end-to-end.
Notebooks 03 and 04 use the same Schrödinger Bridge model; 03 is pedagogical (2D),
04 is the real-data version.

---

## What each notebook produces

### 01 — Paul15 drift field

```
results/01_paul15/
  adata.h5ad
  instability_genes.csv
  regulators.csv
  instability_genes.pdf        ← sensitivity + gene tracks + heatmap
  regulator_summary.pdf        ← 4-panel regulator analysis
  regulator_network.pdf        ← hybrid TF-target graph
```

### 02 — Fourier extension

```
results/02_fourier/
  fourier_model.pt
  power_spectrum.pdf
  real_vs_generated.pdf        ← spectral validation
  band_genes.pdf               ← top genes per frequency band
```

### 03 — Schrödinger Bridge (synthetic)

```
results/03_bridge/
  bridge_model.pt
  bridge_instability.csv       ← forward vs backward per time step
  instability_bridge.pdf       ← instability curves + direction arrows
```

### 04 — Paul15 bridge

```
results/04_paul15_bridge/
  adata_bridge.h5ad
  instability_genes_forward.csv
  instability_genes_backward.csv
  bridge_summary.pdf           ← 7-panel summary
  bridge_genes.pdf             ← gene heatmaps
  regulator_network_forward.pdf
  regulator_network_backward.pdf
```

---

## Data

Paul15 hematopoiesis data (`data/paul15/paul15.h5`) is fetched automatically via
`scanpy.datasets.paul15()` — no manual download needed.

Large generated files (model checkpoints `*.pt`, processed datasets `*.h5ad`) are
excluded from version control. Run the notebooks to regenerate them.
