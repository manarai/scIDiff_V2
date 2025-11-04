
# scQDiff

**scQDiff** learns a **time-conditioned drift field** and extracts **regulatory archetypes** from single-cell data using a Schrödinger-Bridge / score-based formulation—then simulates **forward** (naive→perturbed) and **reverse** transitions. Includes **graph Laplacian** priors, **Fokker–Planck residual**, **AnnData/scVelo/CellRank** support, **fate-conditioned archetypes**, and **trajectory simulation** utilities.

## Install
```bash
# Create conda env
conda create -n scqdiff python=3.10
conda activate scidiff

# Clone the repository
git clone https://github.com/manarai/scQDiff.git
cd sQIDiff

# Install dependencies
pip install -r requirements.txt

# Install in development mode
pip install -e .

```

## AnnData / scVelo / CellRank
- Velocity layers auto-detected: `velocity`, `velocity_s`, `velocity_u` (override with `--vel-layer`).
- Pseudotime keys auto-detected: `rank_pseudotime`, `latent_time`, `dpt_pseudotime` (override with `--ptime-key`).
- Uses `adata.obsp['connectivities']` for Laplacian when present.

**CLI training:**
```bash
python -m scqdiff.pipeline.train_from_anndata
  --h5ad data.h5ad --epochs 200
  --vel-layer velocity
  --ptime-key rank_pseudotime
   --fate-index 0
  --nbins 12
  --rank 3
  --out-prefix results/my_run
  --cpdb-means cpdb/significant_means.parquet
  --cpdb-pvals cpdb/pvalues.parquet
  --time-key bin_or_pseudotime_key
  --comm-knn 20 --comm-gain 0.2

```

## Simulate trajectories
Use `scqdiff.simulate.trajectories` to synthesize forward/reverse dynamics and probe archetype activation:
```python
from scqdiff.simulate.trajectories import euler_integrate, euler_integrate_with_archetype
traj = euler_integrate(model, x0, t0=0.0, t1=1.0, steps=100)      # forward
traj_rev = euler_integrate(model, x1, t0=1.0, t1=0.0, steps=100)  # reverse
traj_mod = euler_integrate_with_archetype(model, x0, patterns, U_t, which=0, lam=0.8)
```
See `examples/05_simulate_trajectories.ipynb`.

## Fate-conditioned archetypes (CellRank)
If your `.h5ad` includes CellRank absorption probabilities (e.g., `obsm['to_fates']`), scQDiff can extract **fate-conditioned** temporal Jacobians. See `examples/04_fate_conditioned_archetypes.ipynb`.
