# scOpAtlas: Stable Operator Atlas

**scOpAtlas** (Stable Operator Atlas) is an analytical layer within the **scjdo** framework that defines cellular states by their **local stability structure** rather than expression patterns alone.

## Overview

Traditional cell atlases define cell states based on gene expression patterns. scOpAtlas adds a **dynamical layer** by characterizing cells according to their **operator regimes**—the local stability properties of the regulatory dynamics governing their behavior.

### Key Concept

> Cell states can be defined by their local operator regime—**stable**, **plastic**, or **unstable**—rather than expression alone.

This reveals:
- **Commitment depth** (how resistant a cell is to perturbation)
- **Plasticity** (how many directions a cell can move)
- **Bifurcation points** (where cell fate decisions occur)
- **Dynamical changes** invisible to expression-based analysis

## Mathematical Foundation

From the learned drift field `f(x,t)` in scjdo, we compute the temporal Jacobian:

```
J(x,t) = ∂f/∂x
```

At each cell location, we extract four key metrics from the eigenvalue spectrum `{λᵢ}` of `J(x,t)`:

### Operator Metrics

| Metric | Formula | Interpretation | Biological Meaning |
|--------|---------|----------------|-------------------|
| **Max Unstable Eigenvalue** | λ_max⁺ = max(Re(λᵢ)) | Detects bifurcations (sign is robust; magnitude/location is not — see gauge note below) | How easily cell escapes current state |
| **Stability Depth** | λ_min⁻ = min(Re(λᵢ))  ·  stability_depth_q = quantile_q(Re(λᵢ)) | Damping strength (extremum drifts with dim — prefer the quantile version for classification and cross-dim comparisons) | Commitment depth, resistance to perturbation |
| **Plasticity Index** | P = #{\|λᵢ\| < ε} / d  (complex modulus) | Fraction of near-neutral modes | Number of accessible directions |
| **Stable Subspace Dimension** | S = #{Re(λᵢ) < 0} | Buffering capacity | Number of stable directions |

> **Gauge freedom.** Snapshot data determines the drift field only up to a
> divergence-free circulation. `stable_dim` and `stability_depth` are
> approximately gauge-invariant; the SIGN of `λ_max⁺` is; its MAGNITUDE at
> bifurcation and the exact spatial location of the crossing are not.
> Pass `symmetric_only=True` to `OperatorMetrics` to eigendecompose
> `(J + Jᵀ)/2` — this restricts to the gauge-robust component and removes
> the modulus/real-part ambiguity, at the cost of dropping the bifurcation
> signal from imaginary parts. See DESIGN.md for the full table.

### Operator Regimes

Four regimes, mutually exclusive by construction. Thresholds are
**data-driven by default** (auto-derived from the dataset's own
spectrum) so results stay comparable across ambient dim.

| Regime | Criteria (defaults, all auto-derived) | Biological Examples |
|--------|---------------------------------------|-------------------|
| **Unstable** | λ_max⁺ > τ_upper (top tercile of λ_max⁺) | Transition states, bifurcations |
| **Plastic** | \|λ_max⁺\| < ε AND plasticity > 0.3 | Progenitor states, decision points |
| **Deeply Stable** | λ_max⁺ ≤ τ_upper AND stability_depth_q < deep_cut | Locked-in fates |
| **Stable** | Everything else | Terminal differentiation, homeostasis |

Pass `thresholds='legacy'` to `atlas.build(...)` to reproduce prior runs
with the old hand-set values (0.1 / 0.05 / -1.0 / 0.3). See DESIGN.md for
the full derivation, including why the extrema-based thresholds collapse
the plastic class at ≥ 50 PCs and how the auto path fixes it.

## Installation

scOpAtlas is included in the scjdo package:

```bash
# Clone repository
git clone https://github.com/manarai/scJDO.git
cd scJDO

# Install dependencies
pip install -r requirements.txt

# Install in development mode
pip install -e .
```

### Dependencies

**Required:**
- `torch` (for autograd Jacobian computation)
- `numpy`
- `anndata`
- `scanpy` (for visualization)
- `matplotlib`, `seaborn`

**Optional:**
- `cellrank` (for trajectory analysis integration)
- `episcanpy` (for ATAC-seq overlay)

## Quick Start

### Python API

```python
import torch
import anndata as ad
from scjdo.models.drift import DriftField
from scjdo.atlas import StableOperatorAtlas

# Load data and trained drift model
adata = ad.read_h5ad("your_data.h5ad")
drift_model = torch.load("your_model.pt")

# Build Stable Operator Atlas
atlas = StableOperatorAtlas(adata, drift_model)
atlas.build()

# Access results
print(atlas.regimes)  # Operator regime labels
print(atlas.metrics)  # Operator metrics

# Validate non-redundancy with cell types
atlas.validate_nonredundancy(celltype_key='cell_type')

# Save results
atlas.save("atlas_results.h5ad")
```

### Command Line Interface

```bash
# Build atlas from trained model
python -m scjdo.atlas.build_atlas_cli \
    --h5ad data.h5ad \
    --model my_model.pt \
    --pseudotime-key pseudotime \
    --celltype-key cell_type \
    --condition-key treatment \
    --out atlas_results.h5ad
```

## Workflow

### 1. Train scjdo Drift Model

First, train a drift field model using scjdo:

```bash
python -m scjdo.pipeline.train_from_anndata \
    --h5ad your_data.h5ad \
    --use-velocity-prior \
    --ptime-key pseudotime \
    --epochs 200 \
    --out-prefix my_model
```

### 2. Build Operator Atlas

```python
from scjdo.atlas import StableOperatorAtlas

atlas = StableOperatorAtlas(
    adata=adata,
    drift_model=drift_model,
    use_rep="X_pca",
    pseudotime_key="pseudotime"
)

atlas.build(
    epsilon=0.1,                    # Threshold for near-neutral modes
    threshold_unstable=0.1,         # Unstable regime threshold
    threshold_plastic=0.05,         # Plastic regime threshold
    threshold_deeply_stable=-1.0,   # Deeply stable threshold
    batch_size=32
)
```

### 3. Validate Non-Redundancy

Critical for demonstrating that operator regimes provide information beyond expression-based cell types:

```python
validation = atlas.validate_nonredundancy(
    celltype_key='cell_type',
    condition_key='condition'
)
```

This tests:
- **Same cell type → different operator regimes** (regime diversity)
- **Same cell type, different conditions → different regime distributions**

### 4. Visualize Results

```python
from scjdo.atlas.visualization import (
    plot_operator_regimes,
    plot_stability_depth_map,
    plot_plasticity_map,
    plot_nonredundancy_comparison
)

# Operator regimes on UMAP
plot_operator_regimes(adata, basis="umap")

# Stability depth map
plot_stability_depth_map(adata, basis="umap")

# Plasticity map
plot_plasticity_map(adata, basis="umap")

# Non-redundancy comparison (critical figure)
plot_nonredundancy_comparison(
    adata,
    celltype_key='cell_type',
    condition_key='condition'
)
```

## Key Features

### 1. Operator Metrics Computation

Compute eigenvalue-derived metrics from the Jacobian of the drift field:

```python
from scjdo.atlas import OperatorMetrics

metrics_computer = OperatorMetrics(drift_model, epsilon=0.1)
metrics = metrics_computer.compute_all_metrics(X, t)

# Access individual metrics
lambda_max = metrics['lambda_max_plus']      # Max unstable eigenvalue
lambda_min = metrics['lambda_min_minus']     # Stability depth
plasticity = metrics['plasticity']           # Plasticity index
stable_dim = metrics['stable_dim']           # Stable subspace dimension
```

### 2. Regime Classification

Classify cells into operator regimes:

```python
from scjdo.atlas import OperatorRegimeClassifier

classifier = OperatorRegimeClassifier(
    threshold_unstable=0.1,
    threshold_plastic=0.05
)

regimes, regime_masks = classifier.classify(metrics)
```

### 3. Condition Comparison

Compare operator regimes across experimental conditions:

```python
comparison = atlas.compare_conditions(
    condition_key='treatment',
    celltype_key='cell_type'
)
```

### 4. Temporal Evolution

Analyze how operator metrics evolve along pseudotime:

```python
from scjdo.atlas.visualization import plot_temporal_evolution

plot_temporal_evolution(
    adata,
    pseudotime_key='pseudotime',
    n_bins=20
)
```

## Biological Applications

### 1. Immune Aging

**Question:** How do immune cells change with age?

**Approach:**
- Compare young vs. old donor samples
- Same cell type (e.g., naïve T cells), different operator regimes
- Hypothesis: Aging deepens stability without major expression changes

```python
atlas.validate_nonredundancy(
    celltype_key='cell_type',
    condition_key='age_group'  # 'young' vs 'old'
)
```

### 2. Drug Resistance

**Question:** What makes cells resistant to treatment?

**Approach:**
- Compare sensitive vs. resistant cell lines
- Identify deeply stable regimes in resistant cells
- Target interventions to shift operator regimes

### 3. Differentiation Commitment

**Question:** When do cells commit to a fate?

**Approach:**
- Track operator metrics along pseudotime
- Identify transition from plastic → stable regimes
- Locate bifurcation points (high λ_max⁺)

### 4. Exhaustion vs. Activation

**Question:** Can exhausted T cells be reactivated?

**Approach:**
- Compare operator regimes of exhausted vs. activated T cells
- Deeply stable exhausted cells → hard to reactivate
- Plastic exhausted cells → reversible

## Output Structure

After building the atlas, results are stored in the AnnData object:

### `adata.obs` (per-cell annotations)

| Column | Description |
|--------|-------------|
| `operator_regime` | Regime label (stable/plastic/unstable/deeply_stable) |
| `lambda_max_plus` | Max Re(λ) per cell |
| `lambda_min_minus` | Min Re(λ) per cell (extremum — dim-drifts; prefer `stability_depth_q`) |
| `stability_depth_q` | Quantile of Re(λ) per cell (dimension-stable stability metric) |
| `plasticity` | Fraction of \|λ\| < ε per cell |
| `stable_dim` | Number of Re(λ) < 0 per cell |
| `boundary_distance` | Clipped, threshold-normalized distance to the classification boundary. **Not a probability** — do not interpret as calibrated confidence. Renamed from `regime_confidence` in Task 2.3; will be replaced by ensemble-concordance confidence when Task 3.2 lands. |

### `adata.uns` (global metadata)

| Key | Description |
|-----|-------------|
| `operator_eigenvalues` | Full eigenvalue spectra (n_cells, n_dims) |

## Visualization Gallery

### 1. Operator Regimes on UMAP

Color cells by operator regime to reveal dynamical structure.

```python
plot_operator_regimes(adata, basis="umap")
```

**Colors:**
- 🟢 Green: Stable
- 🟠 Orange: Plastic
- 🔴 Red: Unstable
- 🔵 Blue: Deeply Stable

### 2. Stability Depth Map

Continuous heatmap showing commitment depth.

```python
plot_stability_depth_map(adata, basis="umap")
```

**Interpretation:**
- Dark regions: Deeply committed states
- Light regions: Shallow commitment, easily perturbed

### 3. Plasticity Map

Highlights decision points and progenitor states.

```python
plot_plasticity_map(adata, basis="umap")
```

**Interpretation:**
- High plasticity: Many accessible directions
- Low plasticity: Constrained motion

### 4. Non-Redundancy Comparison

**Critical figure for publication.** Shows that operator regimes differ across conditions even when expression-based cell types are the same.

```python
plot_nonredundancy_comparison(
    adata,
    celltype_key='cell_type',
    condition_key='condition'
)
```

### 5. Temporal Evolution

Track operator metrics along pseudotime to identify bifurcations and commitment transitions.

```python
plot_temporal_evolution(adata, pseudotime_key='pseudotime')
```

## Validation Strategy

### Non-Redundancy Tests

**Option A (Preferred):** Same cell type → different operator regimes
- Example: Naïve T cells from young vs. old donors
- Same markers, different stability depth

**Option B:** Different cell types → shared operator regime
- Example: Multiple terminal lineages share "maintenance" stability program

### Statistical inference (Task 3.1)

`atlas.validate_nonredundancy(celltype_key, condition_key)` runs the full
inference pipeline in `scjdo/atlas/statistics.py`:

- **Test statistic**: χ² and Cramér's V per cell type on the
  regime × condition contingency table (Cramér's V is a bounded [0, 1]
  effect size — report it alongside the p-value; with 10⁵ cells everything
  is nominally significant).
- **Permutation null**: shuffles CONDITION labels within cell type
  (preserves cell-type composition and regime marginals; destroys the
  condition association). Default 1000 permutations.
- **Depth-matched subsampling**: quantile-bins cells by UMI depth within
  each cell type and equalizes both N and depth distribution per
  condition. Repeat 100× and report the DISTRIBUTION of statistics — not
  a single lucky draw. Regime assignment depends on the drift field,
  which depends on the embedding, which depends on depth — this is the
  step that separates biology from a depth artifact.
- **Multiple testing**: Benjamini–Hochberg q-values across cell types.
- **Bootstrap CIs** on regime fractions per (celltype × condition).
- **Miller-Madow entropy** replaces the plug-in Shannon estimator
  (bias-corrected; small-n effects no longer track sample size).

**Count source** is auto-detected on the AnnData:
`obs['n_counts']` → `obs['total_counts']` → `raw.X` row-sums → `X` row-sums.
Pass `count_key=` to override.

**Acceptance — the negative control** the pipeline must pass before it can
support any headline claim: split a single homogeneous cohort into two
arbitrary halves, label them "A" and "B", run the full pipeline. Median
p across seeded splits must stay well above 0.05. Guarded by
`test_ACCEPTANCE_negative_control` in the test suite.

### Drift-quality gate (Task 3.3)

`atlas.build()` now runs a three-criterion gate before any classification:

1. **Transport error** — pseudotime-holdout SW distance. Bin cells by
   pseudotime; push each bin forward by Δt under the drift; compare the
   pushed distribution to the observed next bin via sliced Wasserstein.
   Pass iff model beats the no-op "stay put" baseline.
2. **Velocity-prior agreement** — mean per-cell cosine similarity
   between drift(x, t) and the RNA-velocity estimate (if present in
   `adata.obsm['velocity']` or `'velocity_pca'`). Not expected to be
   high; expected to be non-random.
3. **Ensemble regime concordance** — reads
   `adata.obs['regime_concordance']` written by `OperatorEnsemble.build`.
   Pass iff the discordant fraction (concordance < 0.5) is under 1/3.

The report is written to `adata.uns['drift_quality']` and travels with
the object. On failure `build()` raises `DriftQualityError` by default;
pass `force=True` to override with a loud banner (the failed report is
still written so downstream code can see the gate was overridden). Turn
the gate off entirely with `quality_check=False`.

Each criterion is skipped gracefully when its inputs aren't available
(no velocity → skip that check; no ensemble → skip that check). The
overall gate passes iff at least one criterion ran and every criterion
that ran passed. On real single-cell data with a legitimate trajectory
you should expect transport to run; on the ideal setup all three run.

### Ensemble across drift-field members (Task 3.2)

`OperatorEnsemble(drift_models)` aggregates regime assignments across
N pre-trained drift fields — the highest-value defense against the
gauge-freedom identifiability problem. Snapshot data determines the
drift only up to a divergence-free (rotational) circulation, so any
single trained field is one arbitrary representative from an
admissible set. An ensemble samples that set.

```python
from scjdo.atlas import OperatorEnsemble

# Load N pre-trained members. Vary at least: seed, network init,
# velocity-prior ON/OFF. Seed-only variation UNDERESTIMATES the
# true gauge spread.
ens = OperatorEnsemble(drift_models=[model_1, model_2, ..., model_N])
ens.build(adata)      # writes adata.obs['ensemble_regime'],
                      # adata.obs['regime_concordance']
ens.print_summary()   # prints modal-regime counts + discordant fractions

# Any single cell's concordance is the fraction of members agreeing with
# its modal assignment — a real per-cell uncertainty quantifier that
# replaces the interim `boundary_distance`.
```

**Report the discordant fraction prominently.** If more than ~1/3 of
cells fall below concordance 0.5, the atlas is not measuring a stable
property — that finding must be stated, not buried.

Once concordance is in `adata.obs['regime_concordance']`,
`atlas.validate_nonredundancy(...)` automatically runs Task 3.1's
inference **twice** — once on all cells, once on concordant-only cells
— and returns both under `statistics` and `statistics_concordant_only`.
Report both.

### Biological Anchoring

Choose ONE strong biological axis:
- **Immune aging** (recommended)
- Exhaustion vs. reversible activation
- Drug resistance vs. sensitivity
- Differentiation commitment

### Chromatin Overlay (Optional)

Use ATAC-seq data to validate:
- Are unstable modes accessible?
- Are stable modes epigenetically reinforced?
- Do resistant states show chromatin locking?

## Advanced Usage

### Custom Thresholds

Adjust classification thresholds based on your data:

```python
atlas.build(
    epsilon=0.15,                    # Wider neutral zone
    threshold_unstable=0.2,          # Higher unstable threshold
    threshold_deeply_stable=-2.0     # Deeper stability requirement
)
```

### Batch Processing

For large datasets, increase batch size:

```python
atlas.build(batch_size=128)  # Process 128 cells at a time
```

### GPU Acceleration

Use GPU for faster computation:

```python
atlas = StableOperatorAtlas(
    adata, drift_model,
    device="cuda"  # Use GPU
)
```

## Troubleshooting

### Issue: "Pseudotime not found"

**Solution:** Compute pseudotime first using Scanpy or CellRank:

```python
import scanpy as sc
sc.tl.diffmap(adata)
sc.tl.dpt(adata)
```

### Issue: "Representation not found"

**Solution:** Compute PCA:

```python
sc.pp.pca(adata, n_comps=50)
```

### Issue: Jacobian computation is slow

**Solutions:**
1. Reduce batch size: `batch_size=16`
2. Use GPU: `device="cuda"`
3. Reduce dimensionality: Use fewer PCs

### Issue: All cells classified as same regime

**Solution:** Adjust thresholds:

```python
atlas.build(
    threshold_unstable=0.05,    # Lower threshold
    threshold_plastic=0.02      # Lower threshold
)
```

## Examples

See `examples/tutorial_scopatlas.py` for a complete tutorial.

## Citation

If you use scOpAtlas in your research, please cite:

> Redd D., Green S., Terooatea T.W. (2026). scJDO: Inferring time-varying dynamical
> operators from single-cell transcriptomic data. *[journal]*.

## Support

For questions and issues:
- GitHub Issues: [github.com/manarai/scJDO/issues](https://github.com/manarai/scJDO/issues)
- Email: tommy.terooatea@byu.edu

## License

MIT License

## Acknowledgments

scOpAtlas builds on the scjdo framework for learning continuous-time cellular dynamics from single-cell data.
