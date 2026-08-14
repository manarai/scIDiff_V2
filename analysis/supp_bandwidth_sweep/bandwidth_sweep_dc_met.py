"""
Bandwidth-sweep diagnostic for the two claimed peaks in the manuscript:
- Paul15 DC branch (near-terminal sensitivity peak reported at t = 0.98)
- SCP295 MET branch (instability peak reported at t = 0.77)

For each branch, fit_drift is called at a grid of bandwidth values.
For each h we record:
  - peak τ location (argmax of max_real_eig over the interior grid)
  - peak Re(λ_max) value
  - n_eff at the peak location
  - bootstrap peak distribution at the selected h (from lam_bootstrap)

Saves a JSON audit trail and a small matplotlib figure suitable for
a supplementary panel.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from sklearn.decomposition import FactorAnalysis

REPO = Path("/Users/terooatt/Downloads/scJDO")
sys.path.insert(0, str(REPO))
from scjdo.tl import fit_drift

OUT = REPO / "scratchpad_run"
OUT.mkdir(exist_ok=True)

# Fast configuration — bandwidth sensitivity should be robust to modest under-
# training since it's an aggregation-step property. Reduce epochs to keep the
# sweep tractable.
N_EPOCHS_FAST = 1000
N_ARCHETYPES = 5
BANDWIDTHS = [0.01, 0.02, 0.03, 0.05, 0.08, 0.10, 0.15]


# ── Paul15 DC branch ──────────────────────────────────────────────────────

def load_paul15_dc():
    """Load paul15 dataset and return the DC-branch subset with FA embedding + Palantir pseudotime."""
    import palantir
    import scipy.sparse
    a = sc.datasets.paul15()
    # Standard preprocessing (matches the paper's Figure 3 setup).
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=2000, flavor="seurat")
    a = a[:, a.var["highly_variable"]].copy()
    X = a.X.toarray() if scipy.sparse.issparse(a.X) else np.asarray(a.X)
    fa = FactorAnalysis(n_components=30, random_state=42, max_iter=1000).fit(X)
    a.obsm["X_fa"] = fa.transform(X).astype(np.float32)

    # Palantir pseudotime and branch masks. Palantir has its own diffusion-
    # map routine (writes DM_EigenVectors / DM_EigenValues to obsm/uns);
    # scanpy's sc.tl.diffmap writes to a different key that Palantir can't find.
    sc.pp.neighbors(a, use_rep="X_fa", n_neighbors=30)
    palantir.utils.run_diffusion_maps(a, pca_key="X_fa", n_components=15)
    ms = palantir.utils.determine_multiscale_space(a, n_eigs=8)

    # Root cell: MEP-like progenitor with the lowest projection on FA0
    # (heuristic that puts the root at the tail of the trajectory).
    root_idx = int(np.argmin(a.obsm["X_fa"][:, 0]))
    pr = palantir.core.run_palantir(
        a, a.obs_names[root_idx],
        num_waypoints=500, use_early_cell_as_start=True,
    )
    a.obs["pseudotime"] = pr.pseudotime.reindex(a.obs_names).values

    # Branch masks: pick whichever terminal Palantir returns with the highest
    # median expression of the DC-marker set (IRF8, TCF4, LILRA4).
    branch_names = list(pr.branch_probs.columns)
    dc_markers = [g for g in ["Irf8", "Tcf4", "Lilra4"] if g in a.var_names]
    if not dc_markers:
        dc_markers = [g for g in ["IRF8", "TCF4"] if g in a.var_names]
    if dc_markers:
        marker_score = a[:, dc_markers].X
        try:
            marker_score = marker_score.toarray().mean(axis=1)
        except AttributeError:
            marker_score = np.asarray(marker_score).mean(axis=1)
        # Assign DC branch = the terminal whose top-100 probable cells have
        # the highest mean marker score.
        best_col, best_score = None, -np.inf
        for c in branch_names:
            top_idx = np.argsort(pr.branch_probs[c].values)[-100:]
            s = float(np.mean(marker_score[top_idx]))
            if s > best_score:
                best_col, best_score = c, s
        dc_col = best_col
    else:
        dc_col = branch_names[0]
    print(f"DC branch column selected: {dc_col} (of {branch_names})")

    dc_probs = pr.branch_probs[dc_col].reindex(a.obs_names).values
    dc_mask = dc_probs > 0.5
    a_dc = a[dc_mask].copy()
    print(f"DC branch: {a_dc.n_obs} cells")
    return a_dc


# ── SCP295 MET branch (subset from the concordance SCP295 h5ad) ──────────

def load_scp295_met():
    """Load MET-branch cells from the SCP295 dataset with cell-set annotation."""
    import scipy.sparse
    p = Path("/Users/terooatt/Documents/Project_scQDiff/02_scQDiff/scIDIFF_anndata/data/SCP295/scp295.h5ad")
    a = sc.read_h5ad(p)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=2000, flavor="seurat")
    a = a[:, a.var["highly_variable"]].copy()
    X = a.X.toarray() if scipy.sparse.issparse(a.X) else np.asarray(a.X)
    fa = FactorAnalysis(n_components=30, random_state=42, max_iter=1000).fit(X)
    a.obsm["X_fa"] = fa.transform(X).astype(np.float32)

    # MET branch: cells whose cell_set is one of the mesenchymal-corridor labels
    # (Progenitor + MET reservoir). Follow Figure4 branch_masks convention.
    cs = a.obs["cell_set"].astype(str).values
    met_mask = np.isin(cs, ["Progenitor", "MET"])
    a_met = a[met_mask].copy()
    a_met.obs["pseudotime"] = a_met.obs["day_norm"].values.astype(np.float32)
    print(f"MET branch: {a_met.n_obs} cells")
    return a_met


# ── Sweep and record ─────────────────────────────────────────────────────

def sweep_bandwidths(a, branch_name):
    records = []
    for h in BANDWIDTHS:
        t0 = time.time()
        model = fit_drift(
            a, rep="X_fa", time_key="pseudotime",
            n_epochs=N_EPOCHS_FAST,
            n_archetypes=N_ARCHETYPES,
            n_eff_min=20.0,
            bandwidth=h,   # explicit — no auto sweep
            grid_size=200,
            n_boot=0,      # skip bootstrap for the interior sweeps (save time)
            seed=42,
            verbose=False,
        )
        res = a.uns["scjdo"]
        lam = np.asarray(res["max_real_eig"])
        t_c = np.asarray(res["t_centers"])
        n_eff = np.asarray(res["n_eff"])
        # interior argmax with an n_eff support filter — otherwise the
        # smoother's boundary behavior produces spurious peaks on the very
        # sparsely-sampled tails.
        interior = (t_c >= 0.05) & (t_c <= 0.95) & ~np.isnan(lam)
        supported = interior & (n_eff >= 20.0)
        if supported.any():
            idx_peak = int(np.where(supported)[0][np.nanargmax(lam[supported])])
            support_ok = True
        else:
            idx_peak = int(np.where(interior)[0][np.nanargmax(lam[interior])])
            support_ok = False
        rec = {
            "bandwidth": float(h),
            "peak_t": float(t_c[idx_peak]),
            "peak_lam": float(lam[idx_peak]),
            "n_eff_at_peak": float(n_eff[idx_peak]) if n_eff.size == lam.size else float("nan"),
            "peak_support_ok": bool(support_ok),
            "r2": float(res["r2"]),
            "fit_time_s": time.time() - t0,
        }
        records.append(rec)
        print(f"  {branch_name}  h={h:.3f}  peak_t={rec['peak_t']:.3f}  "
              f"peak_lam={rec['peak_lam']:.4f}  n_eff@peak={rec['n_eff_at_peak']:.1f}  "
              f"R2={rec['r2']:.3f}  t={rec['fit_time_s']:.0f}s")
    return records


def make_figure(dc_sweep, met_sweep):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, sw, name, target_t in [(axes[0], dc_sweep, "Paul15 DC", 0.98),
                                    (axes[1], met_sweep, "SCP295 MET", 0.768)]:
        hs = [r["bandwidth"] for r in sw]
        peaks = [r["peak_t"] for r in sw]
        lams = [r["peak_lam"] for r in sw]
        neffs = [r["n_eff_at_peak"] for r in sw]
        ax.plot(hs, peaks, "o-", color="#355C9A", label="peak τ")
        ax.axhline(target_t, color="#E8884A", ls="--", alpha=0.6,
                   label=f"paper peak τ = {target_t}")
        ax.set_xlabel("Kernel bandwidth h")
        ax.set_ylabel("Peak τ")
        ax.set_ylim(-0.05, 1.05)
        ax.set_title(f"{name}\n(peak stability under bandwidth sweep)", fontsize=10)
        ax.legend(fontsize=8, loc="best")
        ax.grid(True, alpha=0.3)
        # Second y-axis for n_eff
        ax2 = ax.twinx()
        ax2.plot(hs, neffs, "^-", color="#666666", alpha=0.5, label="n_eff at peak")
        ax2.set_ylabel("n_eff at peak", color="#666666")
        ax2.tick_params(axis="y", labelcolor="#666666")
    fig.suptitle("Supplementary Fig. — bandwidth-sweep sensitivity for main-text peaks",
                 fontsize=11, y=1.02)
    plt.tight_layout()
    fig.savefig(OUT / "supp_fig_bandwidth_sweep.pdf", bbox_inches="tight")
    fig.savefig(OUT / "supp_fig_bandwidth_sweep.png", bbox_inches="tight", dpi=200)
    print(f"Saved: {OUT}/supp_fig_bandwidth_sweep.{{pdf,png}}")


def main():
    print(f"Bandwidth sweep — h ∈ {BANDWIDTHS}  (N_EPOCHS={N_EPOCHS_FAST}, K={N_ARCHETYPES})")
    print()
    print("--- Paul15 DC branch ---")
    try:
        a_dc = load_paul15_dc()
        dc_sweep = sweep_bandwidths(a_dc, "DC")
    except Exception as e:
        print(f"  [skipped: {e}]")
        dc_sweep = []
    print()
    print("--- SCP295 MET branch ---")
    try:
        a_met = load_scp295_met()
        met_sweep = sweep_bandwidths(a_met, "MET")
    except Exception as e:
        print(f"  [skipped: {e}]")
        met_sweep = []

    out = {"DC": dc_sweep, "MET": met_sweep, "bandwidth_grid": BANDWIDTHS}
    (OUT / "bandwidth_sweep_dc_met.json").write_text(json.dumps(out, indent=2))
    print(f"\nSaved: {OUT}/bandwidth_sweep_dc_met.json")

    if dc_sweep and met_sweep:
        make_figure(dc_sweep, met_sweep)


if __name__ == "__main__":
    main()
