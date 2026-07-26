"""Palantir-vs-DPT branch-restricted pseudotime agreement on Paul15.

Reproduces the Fig 3 preprocessing (30-D FA on 1000 HVGs + Palantir), then
computes DPT (k=15 kNN per manuscript Methods) with the same root and reports
the branch-restricted Spearman rank correlation between Palantir and DPT for
Ery, DC, and Mono.

Manuscript claim: Palantir-vs-DPT branch-restricted rank correlation
ρ ≈ 0.08 on the erythroid branch (much lower than the typical Palantir/DPT
agreement of ρ ≳ 0.8 on Paul15 hematopoiesis).
"""
from __future__ import annotations

import json
import os
import warnings
from pathlib import Path

import numpy as np
import scanpy as sc
import palantir
import pandas as pd
from scipy.stats import spearmanr, pearsonr
from sklearn.decomposition import FactorAnalysis
import scipy.sparse as sp

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[3]  # → scJDO repo root
DATA_PATH = REPO / "examples" / "marrow_sample_scseq_counts.h5ad"
OUT_JSON = REPO / "analysis" / "results" / "figure3" / "ery_palantir_dpt_agreement.json"

# Match Figure3_FA.ipynb config
N_HVG = 1000
N_LATENT = 30
FA_MAX_ITER = 5000
N_DM_COMPONENTS = 20
N_WAYPOINTS = 500
SEED = 42

START_CELL = "Run5_164698952452459"
TERMINAL_STATES = pd.Series(
    ["DC", "Mono", "Ery"],
    index=["Run5_131097901611291", "Run5_134936662236454", "Run4_200562869397916"],
)

print(f"Loading {DATA_PATH}")
ad = sc.read(DATA_PATH)
print(f"  {ad.n_obs} cells × {ad.n_vars} genes")

ad.layers["raw_counts"] = ad.X.copy()
sc.pp.normalize_per_cell(ad)
palantir.preprocess.log_transform(ad)
sc.pp.highly_variable_genes(ad, n_top_genes=N_HVG, flavor="cell_ranger")

hvg = ad.var["highly_variable"].values
Xhvg = ad[:, hvg].X
Xhvg = Xhvg.toarray() if sp.issparse(Xhvg) else np.asarray(Xhvg, dtype=np.float32)
fa = FactorAnalysis(n_components=N_LATENT, random_state=SEED, max_iter=FA_MAX_ITER)
ad.obsm["X_fa"] = fa.fit_transform(Xhvg).astype("float32")

# Diffusion maps + kNN graph
palantir.utils.run_diffusion_maps(ad, n_components=N_DM_COMPONENTS, pca_key="X_fa")
palantir.utils.determine_multiscale_space(ad)

# Palantir pseudotime
print("Running Palantir ...")
palantir.core.run_palantir(
    ad, START_CELL, num_waypoints=N_WAYPOINTS, terminal_states=TERMINAL_STATES,
)

# Branch masks
palantir.presults.select_branch_cells(ad, q=0.01, eps=0.01,
                                      masks_key="branch_masks", save_as_df=True)
print("Palantir branch counts:")
for col in ad.obsm["branch_masks"].columns:
    print(f"  {col}: {int(ad.obsm['branch_masks'][col].sum())} cells")

# ---------------------------------------------------------------------------
# DPT computation on the same 30-D FA representation using k=15 kNN
# (matches manuscript Methods description: "diffusion pseudotime (DPT) on a
# k-nearest-neighbor graph (k = 15, Euclidean distance in FA space)")
# ---------------------------------------------------------------------------
print("Building k=15 kNN + DPT root ...")
ad_dpt = ad.copy()
sc.pp.neighbors(ad_dpt, use_rep="X_fa", n_neighbors=15)
sc.tl.diffmap(ad_dpt)
iroot = int(np.where(ad_dpt.obs_names == START_CELL)[0][0])
ad_dpt.uns["iroot"] = iroot
sc.tl.dpt(ad_dpt)

ad.obs["dpt_pseudotime"] = ad_dpt.obs["dpt_pseudotime"].values
print(f"  DPT range: [{ad.obs['dpt_pseudotime'].min():.4f}, "
      f"{ad.obs['dpt_pseudotime'].max():.4f}]")

# ---------------------------------------------------------------------------
# Branch-restricted Spearman ρ between Palantir and DPT
# ---------------------------------------------------------------------------
results = {"root_cell": START_CELL, "n_neighbors": 15, "n_latent": N_LATENT,
           "n_hvg": N_HVG, "seed": SEED}
palantir_pt = ad.obs["palantir_pseudotime"].values
dpt_pt = ad.obs["dpt_pseudotime"].values

whole = spearmanr(palantir_pt, dpt_pt)
results["whole_dataset"] = {
    "n": int(ad.n_obs), "spearman_rho": float(whole.statistic),
    "spearman_p": float(whole.pvalue),
}
print(f"\nWhole-dataset Palantir↔DPT Spearman ρ = {whole.statistic:.3f} "
      f"(n={ad.n_obs})")

print("\n[Variant A] Global DPT restricted to branch mask (matched-root):")
for branch in ["Ery", "DC", "Mono"]:
    mask = ad.obsm["branch_masks"][branch].values.astype(bool)
    n = int(mask.sum())
    if n < 5:
        results[branch + "_global_restricted"] = {"n": n, "note": "too few cells"}
        continue
    rho = spearmanr(palantir_pt[mask], dpt_pt[mask])
    pear = pearsonr(palantir_pt[mask], dpt_pt[mask])
    results[branch + "_global_restricted"] = {
        "n": n,
        "spearman_rho": float(rho.statistic),
        "spearman_p": float(rho.pvalue),
        "pearson_r": float(pear.statistic),
        "pearson_p": float(pear.pvalue),
    }
    print(f"  {branch} (n={n}): ρ = {rho.statistic:.3f} (p={rho.pvalue:.2e})")

print("\n[Variant C] BOTH Palantir and DPT recomputed on the branch subset from "
      "scratch with matched root selection (earliest global-Palantir-pt cell):")
for branch in ["Ery", "DC", "Mono"]:
    mask = ad.obsm["branch_masks"][branch].values.astype(bool)
    n = int(mask.sum())
    if n < 50:
        continue
    sub = ad[mask].copy()
    root_local_name = sub.obs_names[int(np.argmin(sub.obs["palantir_pseudotime"].values))]
    root_local_idx = int(np.where(sub.obs_names == root_local_name)[0][0])
    # Rebuild diffusion structures on the subset
    palantir.utils.run_diffusion_maps(sub, n_components=N_DM_COMPONENTS, pca_key="X_fa")
    palantir.utils.determine_multiscale_space(sub)
    try:
        palantir.core.run_palantir(sub, root_local_name,
                                   num_waypoints=min(N_WAYPOINTS, n // 4))
        palantir_local = sub.obs["palantir_pseudotime"].values
    except Exception as ex:
        print(f"  {branch}: local Palantir failed ({type(ex).__name__}: {ex}); skipping")
        continue
    sc.pp.neighbors(sub, use_rep="X_fa", n_neighbors=15)
    sc.tl.diffmap(sub)
    sub.uns["iroot"] = root_local_idx
    sc.tl.dpt(sub)
    dpt_local = sub.obs["dpt_pseudotime"].values
    finite = np.isfinite(palantir_local) & np.isfinite(dpt_local)
    n_finite = int(finite.sum())
    rho = spearmanr(palantir_local[finite], dpt_local[finite])
    results[branch + "_both_recomputed_on_subset"] = {
        "n": n, "n_finite": n_finite,
        "spearman_rho": float(rho.statistic),
        "spearman_p": float(rho.pvalue),
    }
    print(f"  {branch} (n={n_finite}): ρ = {rho.statistic:.3f} (p={rho.pvalue:.2e})")

print("\n[Variant B] DPT recomputed ON the branch subset with matched root "
      "(earliest Palantir-pseudotime cell in the branch):")
for branch in ["Ery", "DC", "Mono"]:
    mask = ad.obsm["branch_masks"][branch].values.astype(bool)
    n = int(mask.sum())
    if n < 20:
        results[branch + "_branch_subset"] = {"n": n, "note": "too few cells"}
        continue
    sub = ad[mask].copy()
    root_local = int(np.argmin(sub.obs["palantir_pseudotime"].values))
    sc.pp.neighbors(sub, use_rep="X_fa", n_neighbors=15)
    sc.tl.diffmap(sub)
    sub.uns["iroot"] = root_local
    sc.tl.dpt(sub)
    dpt_sub = sub.obs["dpt_pseudotime"].values
    palantir_sub = sub.obs["palantir_pseudotime"].values
    finite = np.isfinite(dpt_sub) & np.isfinite(palantir_sub)
    n_finite = int(finite.sum())
    if n_finite < 20:
        results[branch + "_branch_subset"] = {"n": n, "n_finite": n_finite,
                                              "note": "DPT returned NaN on most cells"}
        continue
    rho = spearmanr(palantir_sub[finite], dpt_sub[finite])
    pear = pearsonr(palantir_sub[finite], dpt_sub[finite])
    results[branch + "_branch_subset"] = {
        "n": n, "n_finite": n_finite,
        "spearman_rho": float(rho.statistic),
        "spearman_p": float(rho.pvalue),
        "pearson_r": float(pear.statistic),
        "pearson_p": float(pear.pvalue),
    }
    print(f"  {branch} (n={n_finite}): ρ = {rho.statistic:.3f} "
          f"(p={rho.pvalue:.2e}); r = {pear.statistic:.3f}")

OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_JSON, "w") as f:
    json.dump(results, f, indent=2)
print(f"\nSaved: {OUT_JSON}")
