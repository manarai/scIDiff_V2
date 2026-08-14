# SCOPAtlas Quick Start Guide

## What is SCOPAtlas?

**SCOPAtlas** (Stable Operator Atlas) defines cellular states by their **local stability structure** rather than expression patterns alone. It provides a dynamical layer of cell identity that is invisible to expression-based atlases.

### Key Difference: scJDO vs SCOPAtlas

| Aspect | scJDO | SCOPAtlas |
|--------|---------|-----------|
| **Purpose** | Learn time-varying drift fields | Define cell states by stability |
| **Jacobian** | Temporal J(t) along trajectories | Per-cell J_ic at each cell |
| **Output** | Regulatory archetypes | Operator regimes & metrics |
| **Analysis** | Temporal evolution | Single-cell classification |
| **Paper focus** | Conserved regulatory grammar | Operator-based state definition |

**Relationship**: SCOPAtlas uses the drift model trained by scJDO to compute per-cell operators.

## Installation

```bash
git clone https://github.com/manarai/scJDO.git
cd scJDO
pip install -e .
```

## Quick Start (5 minutes)

### Step 1: Train scJDO Model

```python
import anndata as ad
from scjdo.models.drift import DriftField
from scjdo.pipeline import train_drift_model

# Load your data
adata = ad.read_h5ad("your_data.h5ad")

# Train drift model
drift_model = train_drift_model(
    adata,
    use_rep='X_pca',
    pseudotime_key='pseudotime',
    n_epochs=100
)

# Save model
drift_model.save("my_drift_model.pt")
```

### Step 2: Build SCOPAtlas

```python
from scjdo.atlas import StableOperatorAtlas

# Initialize atlas
atlas = StableOperatorAtlas(
    adata=adata,
    drift_model=drift_model,
    use_rep="X_pca",
    pseudotime_key="pseudotime"
)

# Build atlas (computes operator metrics)
atlas.build()

# Results are now in adata.obs:
# - 'operator_regime': stable/plastic/unstable/deeply_stable
# - 'lambda_max_plus': Max unstable eigenvalue
# - 'lambda_min_minus': Stability depth
# - 'plasticity': Plasticity index
# - 'stable_dim': Stable subspace dimension
```

### Step 3: Visualize Results

```python
import scanpy as sc

# Compute UMAP
sc.pp.neighbors(adata, use_rep='X_pca')
sc.tl.umap(adata)

# Plot operator regimes
sc.pl.umap(adata, color='operator_regime')

# Plot operator metrics
sc.pl.umap(adata, color=['lambda_max_plus', 'lambda_min_minus', 
                         'plasticity', 'stable_dim'])
```

### Step 4: Operator-Based Clustering

```python
from scjdo.atlas import OperatorClustering

# Initialize clustering
clusterer = OperatorClustering(adata)

# Prepare operator features
clusterer.prepare_operator_features()

# Cluster in operator space
clusterer.cluster_operator_space()

# Or joint clustering (expression + operator)
clusterer.cluster_joint_space(alpha=0.5)

# Compare with expression-based clustering
results = clusterer.compare_clustering_quality(
    methods={
        'Expression': 'leiden',
        'Operator': 'operator_clusters',
        'Joint': 'joint_clusters'
    },
    celltype_key='cell_type'
)
```

### Step 5: Validate Non-Redundancy

```python
# Critical validation: operator regimes ≠ cell types
validation = atlas.validate_nonredundancy(
    celltype_key='cell_type',
    condition_key='condition'  # Optional
)

# This tests:
# 1. Same cell type → different operator regimes
# 2. Different cell types → same operator regime
```

## Complete Workflow

For a complete tutorial with all steps, see:
- **Jupyter notebook**: `examples/scopatlas_complete_workflow.ipynb`
- **Python script**: `examples/tutorial_scopatlas.py`

## Operator Metrics Explained

### Four Key Metrics

| Metric | Symbol | Interpretation | Biological Meaning |
|--------|--------|----------------|-------------------|
| **Max Unstable Eigenvalue** | λ_max⁺ | Bifurcation sensitivity | How easily cell escapes current state |
| **Stability Depth** | λ_min⁻ | Commitment depth | Resistance to perturbation |
| **Plasticity Index** | P | Fraction of neutral modes | Number of accessible directions |
| **Stable Subspace Dim** | S | Buffering capacity | Number of stable directions |

### Four Operator Regimes

| Regime | Criteria | Biological Examples |
|--------|----------|-------------------|
| **Stable** | λ_max⁺ ≤ 0, large S | Terminal differentiation, homeostasis |
| **Plastic** | λ_max⁺ ≈ 0, high P | Progenitor states, decision points |
| **Unstable** | λ_max⁺ > τ | Transition states, bifurcations |
| **Deeply Stable** | Very negative λ_min⁻ | Resistant states, locked-in fates |

## Key Features

### 1. Operator Embedding

Project operators into low-dimensional space:

```python
from scjdo.atlas import compute_operator_embedding

# Compute embedding
embedding = compute_operator_embedding(
    metrics=atlas.metrics,
    method='metrics',  # or 'pca', 'spectrum'
    n_components=4
)

adata.obsm['X_operator'] = embedding
```

### 2. Operator Clustering

Cluster cells by dynamical properties:

```python
from scjdo.atlas import quick_operator_clustering

# Quick clustering
adata = quick_operator_clustering(
    adata,
    method='joint',  # 'operator', 'expression', or 'joint'
    alpha=0.5
)
```

### 3. Non-Redundancy Validation

Demonstrate operator regimes ≠ cell types:

```python
# Validate non-redundancy
validation = atlas.validate_nonredundancy(
    celltype_key='cell_type',
    condition_key='condition'
)

# Visualize
import pandas as pd
import seaborn as sns

crosstab = pd.crosstab(
    adata.obs['cell_type'],
    adata.obs['operator_regime'],
    normalize='index'
)

sns.heatmap(crosstab, annot=True, fmt='.2f')
```

## Common Use Cases

### Use Case 1: Identify Fragile Cell States

```python
# Find unstable cells (sensitive to perturbations)
unstable_mask = adata.obs['operator_regime'] == 'unstable'
unstable_cells = adata[unstable_mask]

print(f"Found {unstable_mask.sum()} unstable cells")
```

### Use Case 2: Compare Aging Effects

```python
# Compare stability across age groups
young_mask = adata.obs['age'] == 'young'
old_mask = adata.obs['age'] == 'old'

young_stability = adata.obs.loc[young_mask, 'lambda_min_minus'].mean()
old_stability = adata.obs.loc[old_mask, 'lambda_min_minus'].mean()

print(f"Young stability: {young_stability:.3f}")
print(f"Old stability: {old_stability:.3f}")
print(f"Difference: {old_stability - young_stability:.3f}")
```

### Use Case 3: Predict Perturbation Response

```python
# Cells with high λ_max⁺ are more sensitive to perturbations
sensitivity = adata.obs['lambda_max_plus']
sensitive_cells = adata[sensitivity > sensitivity.quantile(0.9)]

print(f"Top 10% most sensitive cells: {len(sensitive_cells)}")
```

### Use Case 4: Identify Reprogramming Barriers

```python
# Deeply stable cells are hard to reprogram
deeply_stable_mask = adata.obs['operator_regime'] == 'deeply_stable'
barrier_cells = adata[deeply_stable_mask]

print(f"Cells with reprogramming barriers: {deeply_stable_mask.sum()}")
```

## Troubleshooting

### Issue: "Metric not found in adata.obs"
**Solution**: Run `atlas.build()` first to compute operator metrics.

### Issue: "Representation not found in adata.obsm"
**Solution**: Run `clusterer.prepare_operator_features()` before clustering.

### Issue: Slow Jacobian computation
**Solution**: 
- Reduce batch size: `atlas.build(batch_size=16)`
- Use GPU: `atlas = StableOperatorAtlas(..., device='cuda')`

### Issue: Too many/few clusters
**Solution**: Adjust resolution parameter:
```python
clusterer.cluster_operator_space(resolution=0.5)  # Fewer clusters
clusterer.cluster_operator_space(resolution=2.0)  # More clusters
```

## Advanced Usage

### Custom Operator Thresholds

```python
atlas.build(
    epsilon=0.1,                    # Neutral mode threshold
    threshold_unstable=0.1,         # Unstable regime threshold
    threshold_plastic=0.05,         # Plastic regime threshold
    threshold_deeply_stable=-1.0,   # Deeply stable threshold
    plasticity_threshold=0.3        # Min plasticity for plastic regime
)
```

### Grid Search for Optimal Alpha

```python
alpha_results, best_alpha = clusterer.grid_search_alpha(
    alphas=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    celltype_key='cell_type'
)

print(f"Best alpha: {best_alpha}")
```

### Cluster Purity Analysis

```python
purity = clusterer.compute_cluster_purity(
    cluster_key='operator_clusters',
    celltype_key='cell_type'
)

print(f"Overall purity: {purity['overall_purity']:.3f}")
```

## Citation

If you use SCOPAtlas in your research, please cite:

```bibtex
@article{scJDO2024,
  title={scJDO: Learning Single-Cell Regulatory Dynamics with Hybrid Drift Fields},
  author={Terooatea, Tommy W. and Redd, David},
  journal={bioRxiv},
  year={2024}
}
```

## Support

- **Documentation**: [README.md](README.md)
- **Design document**: [DESIGN.md](DESIGN.md)
- **Examples**: `examples/` directory
- **Issues**: [GitHub Issues](https://github.com/manarai/scJDO/issues)
