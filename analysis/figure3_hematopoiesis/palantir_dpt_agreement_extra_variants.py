"""Extra DPT-variant sweep for Ery branch, hunting for the manuscript's ρ = 0.08.

Prior run (`palantir_dpt_agreement.py`) tested three plausible configurations,
all with matched root + k=15 kNN in FA space, all landed ρ > 0.84 on Ery.
This script tries a broader configuration sweep so the manuscript claim can
either be pinned down or ruled out with wider evidence.

Variants tested (all restricted to the erythroid branch cells):
  D1. DPT on PCA (scanpy default) with k=15 kNN, matched global root
  D2. DPT on PCA with k=15 kNN, matched local root (earliest branch cell)
  D3. DPT on FA, k=30 kNN (scanpy default)
  D4. DPT on FA, k=15 kNN, root = randomly-sampled early-progenitor cell (5x)
  D5. DPT on the Palantir-built diffusion map (reuse dm eigenvectors), matched root
  D6. DPT via palantir directly (`palantir.utils.compute_pseudotime_dpt` if
      exists) as a cross-implementation check
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import scanpy as sc
import palantir
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import FactorAnalysis, PCA
import scipy.sparse as sp

warnings.filterwarnings("ignore")

REPO = Path("/Users/terooatt/Downloads/scJDO")
DATA_PATH = REPO / "examples" / "marrow_sample_scseq_counts.h5ad"
OUT_JSON = REPO / "analysis" / "results" / "figure3" / "ery_palantir_dpt_extra_variants.json"

N_HVG = 1000
N_LATENT = 30
FA_MAX_ITER = 5000
N_DM = 20
N_WAYPOINTS = 500
SEED = 42
START_CELL = "Run5_164698952452459"
TERMINAL_STATES = pd.Series(
    ["DC", "Mono", "Ery"],
    index=["Run5_131097901611291", "Run5_134936662236454", "Run4_200562869397916"],
)

print("Loading …")
ad = sc.read(DATA_PATH)
ad.layers["raw_counts"] = ad.X.copy()
sc.pp.normalize_per_cell(ad)
palantir.preprocess.log_transform(ad)
sc.pp.highly_variable_genes(ad, n_top_genes=N_HVG, flavor="cell_ranger")

hvg = ad.var["highly_variable"].values
Xhvg = ad[:, hvg].X
Xhvg = Xhvg.toarray() if sp.issparse(Xhvg) else np.asarray(Xhvg, dtype=np.float32)
fa = FactorAnalysis(n_components=N_LATENT, random_state=SEED, max_iter=FA_MAX_ITER)
ad.obsm["X_fa"] = fa.fit_transform(Xhvg).astype("float32")
pca = PCA(n_components=N_LATENT, random_state=SEED)
ad.obsm["X_pca"] = pca.fit_transform(Xhvg).astype("float32")

palantir.utils.run_diffusion_maps(ad, n_components=N_DM, pca_key="X_fa")
palantir.utils.determine_multiscale_space(ad)

print("Running Palantir …")
palantir.core.run_palantir(ad, START_CELL, num_waypoints=N_WAYPOINTS,
                           terminal_states=TERMINAL_STATES)
palantir.presults.select_branch_cells(ad, q=0.01, eps=0.01,
                                      masks_key="branch_masks", save_as_df=True)

palantir_pt = ad.obs["palantir_pseudotime"].values
ery_mask = ad.obsm["branch_masks"]["Ery"].values.astype(bool)
iroot_global = int(np.where(ad.obs_names == START_CELL)[0][0])

results = {}


def report(name, rho, note=""):
    results[name] = {"spearman_rho": float(rho), "note": note}
    print(f"  {name:52s} ρ = {rho:+.3f}   {note}")


print("\n--- Extra DPT-variant sweep on the erythroid branch (n=%d) ---" %
      int(ery_mask.sum()))

# D1: DPT on PCA, k=15, global root
sub = ad.copy()
sc.pp.neighbors(sub, use_rep="X_pca", n_neighbors=15)
sc.tl.diffmap(sub)
sub.uns["iroot"] = iroot_global
sc.tl.dpt(sub)
dpt = sub.obs["dpt_pseudotime"].values
rho = spearmanr(palantir_pt[ery_mask], dpt[ery_mask]).statistic
report("D1_DPT-on-PCA-global-root", rho)

# D2: DPT on PCA, k=15, local root (earliest Palantir pt in Ery)
sub = ad[ery_mask].copy()
sc.pp.neighbors(sub, use_rep="X_pca", n_neighbors=15)
sc.tl.diffmap(sub)
sub.uns["iroot"] = int(np.argmin(sub.obs["palantir_pseudotime"].values))
sc.tl.dpt(sub)
dpt = sub.obs["dpt_pseudotime"].values
finite = np.isfinite(dpt) & np.isfinite(sub.obs["palantir_pseudotime"].values)
rho = spearmanr(sub.obs["palantir_pseudotime"].values[finite],
                dpt[finite]).statistic
report("D2_DPT-on-PCA-subset-local-root", rho,
       f"n_finite={int(finite.sum())}")

# D3: DPT on FA, k=30 (scanpy default)
sub = ad.copy()
sc.pp.neighbors(sub, use_rep="X_fa", n_neighbors=30)
sc.tl.diffmap(sub)
sub.uns["iroot"] = iroot_global
sc.tl.dpt(sub)
dpt = sub.obs["dpt_pseudotime"].values
rho = spearmanr(palantir_pt[ery_mask], dpt[ery_mask]).statistic
report("D3_DPT-on-FA-k30-global-root", rho)

# D4: DPT on FA, k=15, root = random early-progenitor cell (5 draws)
early_pool = np.where(palantir_pt < 0.05)[0]  # any near-progenitor cell
rng = np.random.default_rng(SEED)
d4_rhos = []
for i, cand in enumerate(rng.choice(early_pool, size=min(5, len(early_pool)),
                                    replace=False)):
    sub = ad.copy()
    sc.pp.neighbors(sub, use_rep="X_fa", n_neighbors=15)
    sc.tl.diffmap(sub)
    sub.uns["iroot"] = int(cand)
    sc.tl.dpt(sub)
    dpt = sub.obs["dpt_pseudotime"].values
    rho = spearmanr(palantir_pt[ery_mask], dpt[ery_mask]).statistic
    d4_rhos.append(float(rho))
    print(f"  D4_DPT-FA-random-root-{i} (iroot={int(cand)})  ρ = {rho:+.3f}")
results["D4_random-root-samples"] = d4_rhos
results["D4_range"] = [float(min(d4_rhos)), float(max(d4_rhos))]

# D5: DPT reusing Palantir's diffusion map (avoids sc.tl.diffmap rebuild)
try:
    sub = ad.copy()
    sub.uns["iroot"] = iroot_global
    if "DM_EigenVectors_multiscaled" in sub.obsm:
        sub.obsm["X_diffmap"] = sub.obsm["DM_EigenVectors_multiscaled"].astype(
            "float32")
        sc.pp.neighbors(sub, use_rep="X_diffmap", n_neighbors=15)
        sc.tl.diffmap(sub)
        sub.uns["iroot"] = iroot_global
        sc.tl.dpt(sub)
        dpt = sub.obs["dpt_pseudotime"].values
        rho = spearmanr(palantir_pt[ery_mask], dpt[ery_mask]).statistic
        report("D5_DPT-on-palantir-DM-eigenvectors", rho)
    else:
        report("D5_DPT-on-palantir-DM-eigenvectors", float("nan"),
               note="palantir DM eigenvectors not found")
except Exception as e:
    report("D5_DPT-on-palantir-DM-eigenvectors", float("nan"),
           note=f"failed: {type(e).__name__}")

# D6: DPT via scanpy on raw log expression (HVG subset, no dim-reduction)
try:
    sub = ad[:, ad.var["highly_variable"].values].copy()
    sc.pp.neighbors(sub, use_rep="X", n_neighbors=15)
    sc.tl.diffmap(sub)
    sub.uns["iroot"] = iroot_global
    sc.tl.dpt(sub)
    dpt = sub.obs["dpt_pseudotime"].values
    rho = spearmanr(palantir_pt[ery_mask], dpt[ery_mask]).statistic
    report("D6_DPT-on-raw-HVG-log-expr", rho)
except Exception as e:
    report("D6_DPT-on-raw-HVG-log-expr", float("nan"),
           note=f"failed: {type(e).__name__}: {e}")

OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_JSON, "w") as f:
    json.dump(results, f, indent=2)
print(f"\nSaved: {OUT_JSON}")
print("\nAny variant near ρ = 0.08 would justify the manuscript's claim.")
print("Nothing near 0.08 across D1–D6 ⇒ claim not supported.")
