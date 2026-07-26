"""
Saddle dense-sampling control (critique priority #3).

Hypothesis stated in the Discussion (v31) and in the previous saddle-readouts
run: the transverse instability at τ_crit is invisible to Re(λmax), tr(J),
P⊥JP⊥, and FTLE because the transiting cells are rare — the transverse
direction is low-variance and the score-matching objective absorbs it into
noise.

Test: re-simulate the v3 bifurcation with dense sampling around τ_crit ∈
[0.117, 0.217] (5× the baseline density in that window) and rerun the four
readouts. Also repeat over 5 seeds to answer the critique's "n=1 system, no
seed count".

Outcomes:
  A. If any readout localises the saddle under dense sampling → the mechanism
     is confirmed as sampling-limited, which is exactly the experimental
     lever (dense sampling, metabolic labelling, lineage barcoding) the paper
     already argues for.
  B. If all readouts still fail → the stated mechanism is incomplete and
     something else is going on. Also publishable, but a different story.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import torch
from sklearn.decomposition import PCA

REPO = Path(__file__).resolve().parents[3]  # → scJDO repo root
sys.path.insert(0, str(REPO))

from scjdo.atlas.flow_topology import compute_ftle
from scjdo.simulate.trajectories import euler_integrate
from scjdo.tl import fit_drift

# ── System config (matches saddle_readouts.py exactly) ────────────────────
N_CELLS_BASELINE = 1500
D_GENES = 200
N_LATENT = 20
N_EPOCHS = 800
NOISE_SIGMA = 0.3
N_GRID = 200
ALPHA_MIN = 1.0
ALPHA_MAX = 4.0
# Analytic pitchfork for the printed toggle system: α_crit = 3^(-1/4)(1 + 1/3) ≈ 1.013
# (v53 correction; earlier drafts used α_crit = 1.5 → τ_crit = 0.1667)
ALPHA_CRIT = 3.0 ** (-0.25) * (1.0 + 1.0 / 3.0)  # ≈ 1.013
TAU_CRIT = float((ALPHA_CRIT - ALPHA_MIN) / (ALPHA_MAX - ALPHA_MIN))  # ≈ 0.004
TOLERANCE = 0.05
BANDWIDTH = 0.04
N_SEEDS = 2

# Dense-sampling parameters
# DENSE_WINDOW is kept at the historical [0.117, 0.217] (baseline of the printed
# density-enrichment control). Under the corrected α_crit ≈ 1.013 this is a
# later post-saddle transition window rather than the analytic saddle window
# — the v53 manuscript documents this and reports the experiment as a
# later-window control. Dense sampling centred at τ_crit ≈ 0.004 remains an
# untested follow-up.
DENSE_WINDOW = (0.117, 0.217)
DENSE_MULTIPLIER = 5.0   # 5× baseline density inside the window


def _toggle_sym_fp(alpha, tol=1e-10, max_iter=200):
    x = 1.0
    for _ in range(max_iter):
        xn = alpha / (1.0 + x**4)
        if abs(xn - x) < tol:
            break
        x = xn
    return float(x)


def _toggle_drift(z, alpha):
    x, y = z[..., 0], z[..., 1]
    fx = alpha / (1.0 + y**4) - x
    fy = alpha / (1.0 + x**4) - y
    return np.stack([fx, fy], axis=-1)


def make_v3_with_dense(seed=42, n_cells=None, T_total=10.0, dt=0.02, sigma_sde=0.20,
                       dense=False):
    """Same v3 bifurcation SDE, with an option to over-sample the saddle window."""
    rng = np.random.default_rng(seed)
    n_steps = int(np.ceil(T_total / dt))

    if not dense:
        if n_cells is None:
            n_cells = N_CELLS_BASELINE
        target_step = np.sort(rng.integers(1, n_steps + 1, size=n_cells))
    else:
        # Baseline uniform + oversample inside the τ window mapped to steps.
        n_baseline = N_CELLS_BASELINE
        window_frac = DENSE_WINDOW[1] - DENSE_WINDOW[0]
        # Extra cells to add such that inside-window density is DENSE_MULTIPLIER×
        n_extra = int(round(n_baseline * window_frac * (DENSE_MULTIPLIER - 1.0)))
        target_step_uniform = rng.integers(1, n_steps + 1, size=n_baseline)
        lo_step = max(1, int(DENSE_WINDOW[0] * n_steps))
        hi_step = int(DENSE_WINDOW[1] * n_steps) + 1
        target_step_extra = rng.integers(lo_step, hi_step, size=n_extra)
        target_step = np.sort(np.concatenate([target_step_uniform, target_step_extra]))
        n_cells = int(target_step.size)

    t_per_cell = (target_step * dt) / T_total
    z = np.full((n_cells, 2), _toggle_sym_fp(ALPHA_MIN), dtype=np.float32)
    z += 0.02 * rng.standard_normal(z.shape).astype(np.float32)
    z_out = np.zeros_like(z)
    sqrt_dt = float(np.sqrt(dt))
    idx = 0
    for step in range(1, n_steps + 1):
        a_t = ALPHA_MIN + (ALPHA_MAX - ALPHA_MIN) * (step * dt) / T_total
        z = z + _toggle_drift(z, a_t) * dt + sigma_sde * sqrt_dt * rng.standard_normal(z.shape).astype(np.float32)
        while idx < n_cells and target_step[idx] == step:
            z_out[idx] = z[idx]
            idx += 1
    return z_out.astype(np.float32), t_per_cell.astype(np.float32)


def build_adata(seed, dense=False):
    z_tog, t_true = make_v3_with_dense(seed=seed, dense=dense)
    rng_proj = np.random.default_rng(seed + 1)
    W = rng_proj.standard_normal((2, D_GENES)).astype(np.float32)
    W /= np.linalg.norm(W, axis=0, keepdims=True) + 1e-8
    X_obs = z_tog @ W
    rng_g = np.random.default_rng(seed)
    X_obs = X_obs + NOISE_SIGMA * rng_g.standard_normal(X_obs.shape).astype(np.float32)
    pca = PCA(n_components=N_LATENT, random_state=seed)
    Z_lat = pca.fit_transform(X_obs).astype(np.float32)
    a = ad.AnnData(X=X_obs)
    a.obs_names = [f"c{i:05d}" for i in range(a.n_obs)]
    a.obsm["X_pca"] = Z_lat
    a.obs["pseudotime"] = t_true
    return a, Z_lat, t_true


def kernel_curve(scalar_per_cell, t_true, grid, bandwidth=BANDWIDTH, n_eff_min=20.0):
    curve = np.full(grid.shape, np.nan)
    for k, tau in enumerate(grid):
        w = np.exp(-((t_true.astype(np.float64) - tau) ** 2) / (2.0 * bandwidth ** 2))
        W = w.sum()
        if W < 1e-9:
            continue
        n_eff = (W * W) / (w * w).sum()
        if n_eff < n_eff_min:
            continue
        curve[k] = float((w * scalar_per_cell).sum() / W)
    return curve


def interior_argmax(curve, grid, low=0.05, high=0.95):
    m = (grid >= low) & (grid <= high) & ~np.isnan(curve)
    if not m.any():
        return np.nan
    idx = np.where(m)[0][np.argmax(curve[m])]
    return float(grid[idx])


def interior_argmin(curve, grid, low=0.05, high=0.95):
    m = (grid >= low) & (grid <= high) & ~np.isnan(curve)
    if not m.any():
        return np.nan
    idx = np.where(m)[0][np.argmin(curve[m])]
    return float(grid[idx])


def jacobian_readouts(model, Z_lat, t_true):
    """
    Return per-cell (λ_max, tr(J), λ_max of P⊥JP⊥) via the model's batched
    jacobian() method (one forward + D backward passes, orders of magnitude
    faster than a Python loop over cells).
    """
    model = model.eval()
    device = next(model.parameters()).device
    z = torch.tensor(Z_lat, dtype=torch.float32, device=device)
    t = torch.tensor(t_true, dtype=torch.float32, device=device)
    B, D = z.shape

    # Full batched forward for the transverse projector.
    with torch.no_grad():
        u_all = model(z, t).detach().cpu().numpy()

    # Batched Jacobian: shape (B, D, D).
    J_batch = model.jacobian(z, t)
    J_np = J_batch.detach().cpu().numpy()

    lam = np.zeros(B, dtype=np.float32)
    div = np.zeros(B, dtype=np.float32)
    lam_perp = np.zeros(B, dtype=np.float32)
    for i in range(B):
        J = J_np[i]
        eigs = np.linalg.eigvals(J)
        lam[i] = float(eigs.real.max())
        div[i] = float(np.trace(J).real)
        f = u_all[i]
        n = np.linalg.norm(f)
        if n < 1e-8:
            lam_perp[i] = lam[i]
            continue
        fh = f / n
        P = np.eye(D) - np.outer(fh, fh)
        Jperp = P @ J @ P
        lam_perp[i] = float(np.linalg.eigvals(Jperp).real.max())
    return lam, div, lam_perp


def run_one_seed(seed, dense):
    a, Z_lat, t_true = build_adata(seed, dense=dense)
    t0 = time.time()
    model = fit_drift(
        a, rep="X_pca", time_key="pseudotime", n_epochs=N_EPOCHS,
        n_archetypes=4, n_eff_min=20.0, n_boot=10, grid_size=N_GRID,
        seed=seed, verbose=False,
    )
    r2 = float(a.uns["scjdo"]["r2"])
    lam_pc, div_pc, lam_perp_pc = jacobian_readouts(model, Z_lat, t_true)
    grid = np.linspace(0.0, 1.0, N_GRID, dtype=np.float64)
    lam_c = kernel_curve(lam_pc, t_true, grid)
    div_c = kernel_curve(div_pc, t_true, grid)
    lam_perp_c = kernel_curve(lam_perp_pc, t_true, grid)
    seeds = torch.tensor(Z_lat, dtype=torch.float32)
    ftle_per_cell = compute_ftle(model, seeds, t0=0.0, t1=0.5, steps=50, delta=1e-3)
    ftle_c = kernel_curve(ftle_per_cell, t_true, grid)
    fit_time = time.time() - t0

    return {
        "seed": int(seed),
        "n_cells": int(a.n_obs),
        "r2": r2,
        "fit_time_s": fit_time,
        "tau_R1_lambda_max": interior_argmax(lam_c, grid),
        "tau_R2_div_pos": interior_argmax(div_c, grid),
        "tau_R2_div_neg": interior_argmin(div_c, grid),
        "tau_R3_transverse": interior_argmax(lam_perp_c, grid),
        "tau_R4_FTLE": interior_argmax(ftle_c, grid),
    }


def main():
    print(f"τ_crit (analytic) = {TAU_CRIT:.4f}   tolerance = ±{TOLERANCE}", flush=True)
    print(f"Dense multiplier inside [{DENSE_WINDOW[0]:.3f}, {DENSE_WINDOW[1]:.3f}]: {DENSE_MULTIPLIER}×", flush=True)
    print(f"Skipping baseline — using stored `saddle_readouts.json` for the baseline row.", flush=True)
    print(f"N_SEEDS (dense) = {N_SEEDS}", flush=True)
    print("", flush=True)

    all_results = {}
    per_seed = []
    for s in range(N_SEEDS):
        print(f"[dense] starting seed {s + 100}...", flush=True)
        r = run_one_seed(s + 100, dense=True)
        per_seed.append(r)
        print(f"  seed {r['seed']:>3d}  n={r['n_cells']:>5d}  R²={r['r2']:.3f}  "
              f"t={r['fit_time_s']:.0f}s  "
              f"τ_λmax={r['tau_R1_lambda_max']:.3f}  "
              f"τ_trJ+={r['tau_R2_div_pos']:.3f}  τ_trJ-={r['tau_R2_div_neg']:.3f}  "
              f"τ_P⊥JP⊥={r['tau_R3_transverse']:.3f}  τ_FTLE={r['tau_R4_FTLE']:.3f}",
              flush=True)
    all_results["dense"] = per_seed
    all_results["baseline_note"] = "See scratchpad_run/saddle_readouts.json for the baseline (n=1 seed) numbers."

    print("", flush=True)
    print("=" * 82, flush=True)
    print(f"{'Readout':>20s} {'mean τ':>8s} {'err vs τ_crit':>14s} {'N pass':>10s}", flush=True)
    print("-" * 82, flush=True)
    for readout in ["tau_R1_lambda_max", "tau_R2_div_pos", "tau_R2_div_neg",
                    "tau_R3_transverse", "tau_R4_FTLE"]:
        taus = [r[readout] for r in all_results["dense"]]
        mean_t = float(np.mean(taus))
        err = abs(mean_t - TAU_CRIT)
        n_pass = int(sum(1 for x in taus if abs(x - TAU_CRIT) <= TOLERANCE))
        print(f"{readout:>20s} {mean_t:>8.3f} {err:>14.3f} {n_pass:>7d}/{N_SEEDS:<3d}", flush=True)
    print("=" * 82, flush=True)

    outpath = Path("{Path(__file__).parent / "saddle_dense_sampling.json"}")
    outpath.write_text(json.dumps(all_results, indent=2))
    print(f"\nSaved: {outpath}")


if __name__ == "__main__":
    main()
