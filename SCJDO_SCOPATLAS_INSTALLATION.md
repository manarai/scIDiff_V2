# Installing and Configuring scJDO with scOpAtlas

This guide installs the unified **scJDO** package and configures the integrated **scOpAtlas** module. scJDO is the dynamics engine: it learns a smooth drift field in a latent cell-state representation. scOpAtlas is the interpretation layer: it computes local Jacobian-derived metrics, quality-gates the inferred dynamics, and writes operator-state annotations back to the same `AnnData` object.

> **Scope and interpretation.** The package can be installed and tested on any supported environment. Meaningful biological interpretation requires a well-sampled system with a defensible trajectory, perturbation, spatial-gradient, or temporal signal and a drift model that passes its quality checks. Do not interpret an operator regime as a biological fact solely because the code executed.

## 1. Requirements

| Component | Minimum | Recommended for real analyses | Purpose |
|---|---:|---:|---|
| Python | 3.9 | 3.11 | Package runtime |
| RAM | 16 GB | 32–64 GB | AnnData, neighborhood graph, Jacobian batches |
| CPU | 4 cores | 8+ cores | Preprocessing, drift fitting, operator metrics |
| GPU | Optional | CUDA-capable GPU | Faster neural drift-model training |
| Input | AnnData with counts | AnnData with count layer, covariates, and known/estimated root | Trajectory preparation and modeling |

The verified package declares the required Scanpy/AnnData, PyTorch, scVelo, CellRank, `igraph`, and `leidenalg` dependencies in `pyproject.toml`. `igraph` and `leidenalg` are required for the bundled Leiden-based operator clustering workflow.

## 2. Create an isolated environment

Conda is recommended for reproducible scientific environments. Create the environment, activate it, and update build tooling:

```bash
conda create -n scjdo python=3.11 -y
conda activate scjdo
python -m pip install --upgrade pip setuptools wheel
```

For a CUDA-enabled PyTorch installation, install the PyTorch build appropriate to the target CUDA runtime **before** installing scJDO. Follow the PyTorch installation selector for the exact command matching the target system.

## 3. Install the verified source tree

Clone the repository and install the package in editable mode with development dependencies:

```bash
git clone https://github.com/manarai/scJDO.git
cd scJDO
python -m pip install -e '.[dev]'
```

If you are using the verified update prepared during the audit, copy its modified files into a fresh clone or use the supplied update directory, then install from that directory.

Confirm that the package and both command-line entry points resolve:

```bash
python -c "import scjdo; print(scjdo.__version__ if hasattr(scjdo, '__version__') else 'scjdo import OK')"
scjdo --help
scjdo-atlas --help
pytest -q
```

The verified test baseline is **153 passed and 7 skipped** for the full suite, including **50 scOpAtlas-specific tests**. A different result may indicate platform-specific dependency incompatibility or local source changes.

## 4. Install optional analysis capabilities

The core PCA/DPT workflow has no optional requirement beyond the base package. Install optional packages only when the corresponding workflow is needed.

| Use case | Installation command | Configuration |
|---|---|---|
| Batch-aware latent representation | `python -m pip install scvi-tools` | Use `latent='scvi'` and pass `batch_key` to `sjd.pp.prepare_trajectory`. |
| Branch-aware pseudotime | `python -m pip install palantir` | Use `pseudotime_method='palantir'`; pass resulting branch probabilities to downstream drift fitting when appropriate. |
| Enrichr/GSEA pathway analysis | `python -m pip install gseapy` | Use after exporting driver genes from scJDO/scOpAtlas. |
| Notebook authoring | `python -m pip install jupyterlab ipykernel` | Register the environment as a Jupyter kernel if desired. |

Registering the activated environment in Jupyter is optional but convenient:

```bash
python -m ipykernel install --user --name scjdo --display-name "Python (scJDO)"
```

## 5. Prepare an AnnData object

Start with raw counts when possible. The standard PCA/DPT pathway normalizes data, selects highly variable genes, computes PCA and neighbors, optionally creates UMAP, and stores trajectory metadata in `adata.uns['scjdo_prep']`.

```python
import scanpy as sc
import scjdo as sjd

adata = sc.read_h5ad("input_counts.h5ad")

# Example root choice: a known stem/progenitor/early-condition label.
sjd.pp.prepare_trajectory(
    adata,
    groupby="cell_type",
    root="HSC",
    latent="pca",
    n_hvg=2_000,
    n_pcs=50,
    n_neighbors=15,
    pseudotime_method="dpt",
    time_key="pseudotime",
)
```

After successful preparation, verify the expected artifacts before fitting dynamics:

```python
assert "X_pca" in adata.obsm
assert "pseudotime" in adata.obs
assert adata.uns["scjdo_prep"]["rep"] == "X_pca"
```

### Batch-aware configuration

If batches dominate the PCA structure, use a batch-aware latent space rather than trying to correct the inferred Jacobians after the fact. The following uses scVI and preserves the declared batch covariate:

```python
sjd.pp.prepare_trajectory(
    adata,
    groupby="cell_type",
    root="HSC",
    latent="scvi",
    batch_key="donor",
    n_latent=20,
    n_scvi_epochs=400,
    pseudotime_method="dpt",
)

representation_key = "X_scvi"
```

A batch-aware representation mitigates technical separation but does **not** prove that residual flow is biological. Inspect batch mixing, label preservation, and whether pseudotime is coherent before drift fitting.

## 6. Fit the scJDO dynamics engine

Fit a drift model only after the trajectory and representation are inspected. Use the same `time_key` and representation that were prepared above.

```python
model = sjd.tl.fit_drift(
    adata,
    time_key="pseudotime",
    use_rep=adata.uns["scjdo_prep"]["rep"],
    n_epochs=1_000,
    n_archetypes=4,
    verbose=True,
)
```

For a fast software smoke test, use a small epoch count. For biological analysis, increase the training budget and retain model diagnostics. A run that produces a file is not necessarily converged or biologically identifiable.

## 7. Build the scOpAtlas layer

`StableOperatorAtlas` consumes the trained drift model and writes results directly to `adata`.

```python
from scjdo.atlas import StableOperatorAtlas

atlas = StableOperatorAtlas(
    adata=adata,
    drift_model=model,
    use_rep=adata.uns["scjdo_prep"]["rep"],
    pseudotime_key="pseudotime",
)

atlas.build(
    epsilon=0.1,
    batch_size=32,
    compute_confidence=True,
    quality_check=True,
    force=False,
)
```

The default `quality_check=True` is intentional. If the learned field fails the drift-quality gate, `atlas.build()` raises an error while writing the diagnostic report to `adata.uns['drift_quality']`. Use `force=True` only for exploration after recording why the gate failed; never treat a forced result as validated evidence.

### Main scOpAtlas outputs

| Location | Typical fields | Interpretation |
|---|---|---|
| `adata.obs` | `lambda_max_plus`, `stability_depth_q`, `plasticity`, `operator_regime`, `boundary_distance` | Cell-level operator metrics and regime labels |
| `adata.obsm` | `X_operator` after feature preparation; `X_joint` after joint clustering | Operator-only and joint feature spaces |
| `adata.uns` | `drift_quality`, optional eigenvalue spectra, atlas metadata | Quality and provenance information |

A categorical regime is a concise annotation, not a substitute for examining the continuous metrics, confidence, quality report, and biological context.

## 8. Cluster and visualize operator outputs

```python
from scjdo.atlas.clustering import OperatorClustering

clusterer = OperatorClustering(adata)
clusterer.prepare_operator_features()
clusterer.cluster_operator_space(resolution=1.0, key_added="operator_clusters")
clusterer.cluster_joint_space(
    expression_rep=adata.uns["scjdo_prep"]["rep"],
    alpha=0.5,
    resolution=1.0,
    key_added="joint_clusters",
)

sc.pl.umap(
    adata,
    color=["cell_type", "operator_regime", "plasticity", "joint_clusters"],
    ncols=2,
)
```

The joint representation concatenates standardized expression and operator blocks. `alpha` specifies the relative operator contribution to squared distance; values outside `[0, 1]` are rejected.

## 9. Save results and preserve provenance

Store the enriched AnnData object and a separate model checkpoint. Preserve the software version, parameter dictionary, representation key, pseudotime key, and random seed with each analysis.

```python
adata.write_h5ad("results/scjdo_scopatlas.h5ad")

# Save the model only if the checkpoint is trusted when reloaded.
import torch
torch.save({"state_dict": model.state_dict(), "cfg": model.cfg}, "results/drift_model.pt")
```

When loading a local, trusted checkpoint produced by the same workflow under current PyTorch releases, use the explicit trusted-object setting:

```python
checkpoint = torch.load("results/drift_model.pt", map_location="cpu", weights_only=False)
```

Never set `weights_only=False` when loading an untrusted checkpoint from an unknown source, because PyTorch pickle deserialization can execute arbitrary code.

## 10. Troubleshooting

| Symptom | Likely cause | Recommended action |
|---|---|---|
| No pseudotime column | `root` and `groupby` were not both supplied or labels do not match | Confirm the root label and run trajectory preparation again. |
| Quality gate fails | Weak trajectory, batch-driven representation, undertrained model, or unsupported velocity/transport assumptions | Inspect `adata.uns['drift_quality']`; improve representation/trajectory and retrain rather than immediately using `force=True`. |
| All cells have one regime | Undertrained model, overly coarse data, static biology, inappropriate thresholds, or genuinely homogeneous dynamics | First increase training/validate field; then inspect continuous metrics and thresholds. |
| Leiden import/cluster failure | Missing graph dependencies | Reinstall the current package metadata; confirm `igraph` and `leidenalg` import. |
| H5AD save error | Unsupported custom Python object in `adata.uns` | Use the updated repository code and store serializable arrays/tables only. |
| No instability-driver genes | No region passes the sensitivity threshold | Report the absence; do not force a gene list. Consider whether the system contains a transition and whether the field passed quality checks. |

## 11. Minimal reproducibility checklist

Before interpreting or sharing an atlas, record the following:

1. Dataset accession/source, preprocessing choices, and genes retained.
2. Representation backend and batch-correction configuration.
3. Root definition, pseudotime method, and relevant trajectory diagnostics.
4. Drift-model architecture, epochs, seed, and fit diagnostics.
5. scOpAtlas thresholds, quality-gate result, and whether `force=True` was used.
6. Regime distributions plus continuous operator metrics, not only UMAP colors.
7. Exact package version/commit and exported `AnnData` file.

## 12. Useful commands

```bash
# Run the complete tests
pytest -q

# Run only scOpAtlas regression tests
pytest tests/test_scopatlas.py -q

# Verify both CLIs
scjdo --help
scjdo-atlas --help

# Execute a tutorial headlessly
MPLBACKEND=Agg jupyter nbconvert --to notebook --execute \
  examples/05_scopatlas_complete_workflow.ipynb
```

## References

[1] [scJDO source repository](https://github.com/manarai/scJDO)

[2] [PyTorch documentation: serialization semantics](https://pytorch.org/docs/stable/notes/serialization.html)

[3] [scVI-tools documentation](https://docs.scvi-tools.org/)

[4] [Palantir documentation](https://palantir.readthedocs.io/)
