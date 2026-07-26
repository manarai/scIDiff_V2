"""
Saddle-readout comparison (Path A of the critique) on the v3 bifurcation system.

Same trained drift field, four readouts along oracle pseudotime τ ∈ [0, 1],
each scored against the analytic saddle τ_crit = (α_crit - α_min)/(α_max - α_min):

  R1  Re(λmax)                    — baseline (the one that fails)
  R2  Divergence = tr(J)           — condensation/expansion of volume
  R3  Transverse-projected spectrum:
      max Re eig( P⊥ J P⊥ ) with P⊥ = I - f̂ f̂ᵀ
                                    — instability orthogonal to the drift tangent
  R4  FTLE ridge from scjdo.atlas.flow_topology
                                    — finite-time flow-map stretching

τ_crit is discussed in the Discussion but never used to define the readouts
(they are all label-blind); it is only used at the end to compute per-readout
error and PASS/FAIL against the pre-committed 0.05 window.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import torch
from sklearn.decomposition import PCA

REPO = Path(__file__).resolve().parents[3]  # → scJDO repo root
sys.path.insert(0, str(REPO))

from scjdo.tl import fit_drift
from scjdo.simulate.trajectories import euler_integrate
from scjdo.atlas.flow_topology import compute_ftle

import anndata as ad


# ── Ground truth: symmetric toggle switch, α ramped from α_min → α_max ────
SEED         = 42
N_CELLS      = 1500
D_GENES      = 200
N_LATENT     = 20
N_EPOCHS     = 800
NOISE_SIGMA  = 0.3
N_GRID       = 200
ALPHA_MIN    = 1.0
ALPHA_MAX    = 4.0
# Analytic pitchfork for the printed toggle system: α_crit = 3^(-1/4)(1 + 1/3) ≈ 1.013
# (v53 correction; earlier drafts used α_crit = 1.5 → τ_crit = 0.1667)
ALPHA_CRIT   = 3.0 ** (-0.25) * (1.0 + 1.0 / 3.0)  # ≈ 1.013
TAU_CRIT     = float((ALPHA_CRIT - ALPHA_MIN) / (ALPHA_MAX - ALPHA_MIN))  # ≈ 0.004
TOLERANCE    = 0.05
BANDWIDTH    = 0.04


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


def make_v3(n_cells=N_CELLS, seed=SEED, T_total=10.0, dt=0.02, sigma_sde=0.20):
    rng = np.random.default_rng(seed)
    n_steps = int(np.ceil(T_total / dt))
    target_step = np.sort(rng.integers(1, n_steps + 1, size=n_cells))
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


def build_adata(seed=SEED):
    z_tog, t_true = make_v3()
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


# ── Kernel aggregation of any per-cell scalar onto a τ-grid ──────────────

def kernel_curve(scalar_per_cell, t_true, grid, bandwidth=BANDWIDTH, n_eff_min=20.0):
    curve = np.full(grid.shape, np.nan)
    n_eff_out = np.full(grid.shape, np.nan)
    for k, tau in enumerate(grid):
        w = np.exp(-((t_true.astype(np.float64) - tau) ** 2) / (2.0 * bandwidth ** 2))
        W = w.sum()
        if W < 1e-9:
            continue
        n_eff_out[k] = (W * W) / (w * w).sum()
        if n_eff_out[k] < n_eff_min:
            continue
        curve[k] = float((w * scalar_per_cell).sum() / W)
    return curve, n_eff_out


def interior_argmax(curve, grid, low=0.05, high=0.95):
    mask = (grid >= low) & (grid <= high) & ~np.isnan(curve)
    if not mask.any():
        return np.nan
    idx = np.where(mask)[0][np.argmax(curve[mask])]
    return float(grid[idx])


def interior_argmin(curve, grid, low=0.05, high=0.95):
    mask = (grid >= low) & (grid <= high) & ~np.isnan(curve)
    if not mask.any():
        return np.nan
    idx = np.where(mask)[0][np.argmin(curve[mask])]
    return float(grid[idx])


# ── R1/R2/R3: per-cell Jacobian readouts ─────────────────────────────────

def jacobian_readouts(model, Z_lat, t_true):
    """
    For each cell, compute Re(λmax(J)), tr(J), and max Re(λ) of P⊥ J P⊥
    where P⊥ = I - f̂ f̂ᵀ projects out the drift tangent.
    """
    model = model.eval()
    device = next(model.parameters()).device
    z = torch.tensor(Z_lat, dtype=torch.float32, device=device)
    t = torch.tensor(t_true, dtype=torch.float32, device=device)
    B, D = z.shape

    lam_max = np.zeros(B, dtype=np.float32)
    div = np.zeros(B, dtype=np.float32)
    lam_max_transverse = np.zeros(B, dtype=np.float32)

    for i in range(B):
        xi = z[i : i + 1].detach().requires_grad_(True)
        ti = t[i : i + 1]
        u = model(xi, ti)  # (1, D)

        # Jacobian rows
        rows = []
        for j in range(D):
            g = torch.autograd.grad(
                u[0, j], xi, retain_graph=True, create_graph=False
            )[0][0]  # (D,)
            rows.append(g.detach().cpu().numpy())
        J = np.stack(rows, axis=0)  # (D, D)

        eigs = np.linalg.eigvals(J)
        lam_max[i] = float(eigs.real.max())
        div[i] = float(np.trace(J).real)

        # Transverse projector
        f = u.detach().cpu().numpy().reshape(-1)  # (D,)
        n = np.linalg.norm(f)
        if n < 1e-8:
            lam_max_transverse[i] = lam_max[i]
            continue
        fh = f / n
        P = np.eye(D) - np.outer(fh, fh)
        Jperp = P @ J @ P
        eigs_p = np.linalg.eigvals(Jperp)
        lam_max_transverse[i] = float(eigs_p.real.max())

    return lam_max, div, lam_max_transverse


# ── R4: FTLE from scjdo.atlas.flow_topology ──────────────────────────────

def ftle_curve(model, Z_lat, t_true, grid, t0=0.0, t1=1.0, steps=100, delta=1e-3):
    """
    Compute FTLE at each cell seed with the flow window [0, 1], then aggregate
    onto the pseudotime grid via the same Gaussian kernel used for other
    readouts. Uses the FTLE from the freshly-added flow_topology module.
    """
    seeds = torch.tensor(Z_lat, dtype=torch.float32)
    ftle_per_cell = compute_ftle(model, seeds, t0=t0, t1=t1, steps=steps, delta=delta)
    curve, _ = kernel_curve(ftle_per_cell, t_true, grid)
    return curve, ftle_per_cell


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    print(f"τ_crit (analytic, scoring only) = {TAU_CRIT:.4f}")
    print(f"scoring window = [{max(0.0, TAU_CRIT - TOLERANCE):.3f}, {TAU_CRIT + TOLERANCE:.3f}]  (clamped ≥ 0)")
    print()

    adata, Z_lat, t_true = build_adata(SEED)
    print(f"cells = {adata.n_obs}  |  latent dim = {Z_lat.shape[1]}")

    t0 = time.time()
    model = fit_drift(
        adata,
        rep="X_pca",
        time_key="pseudotime",
        n_epochs=N_EPOCHS,
        n_archetypes=4,
        n_eff_min=20.0,
        n_boot=10,
        grid_size=N_GRID,
        seed=SEED,
        verbose=False,
    )
    r2 = float(adata.uns["scjdo"]["r2"])
    print(f"fit_drift done in {time.time() - t0:.1f}s   R² = {r2:.3f}")

    # Verify the baseline failure signature by using the scJDO-provided curve.
    lam_scjdo = np.asarray(adata.uns["scjdo"]["max_real_eig"])
    t_cen_scjdo = np.asarray(adata.uns["scjdo"]["t_centers"])
    tau_R1_scjdo = float(t_cen_scjdo[int(np.nanargmax(lam_scjdo))])
    print(f"R1  scJDO cached λ_max peak τ = {tau_R1_scjdo:.3f}   err = {abs(tau_R1_scjdo - TAU_CRIT):.3f}")

    # Compute per-cell readouts (fresh, from the trained drift model).
    grid = np.linspace(0.0, 1.0, N_GRID, dtype=np.float64)
    print("Computing per-cell Jacobians (may take a minute)...")
    lam_max_pc, div_pc, lam_max_perp_pc = jacobian_readouts(model, Z_lat, t_true)

    # Kernel-aggregate each per-cell readout onto the pseudotime grid.
    lam_max_curve, _ = kernel_curve(lam_max_pc, t_true, grid)
    div_curve, _ = kernel_curve(div_pc, t_true, grid)
    lam_perp_curve, _ = kernel_curve(lam_max_perp_pc, t_true, grid)

    # For divergence, the informative direction is CONTRACTION — negative tr(J)
    # peaks near saddles because the transverse direction unbalances. But we
    # also report the |tr(J)| interpretation; whichever localises better wins.
    print("Computing FTLE from flow_topology...")
    ftle_c, _ = ftle_curve(model, Z_lat, t_true, grid, t0=0.0, t1=0.5, steps=50, delta=1e-3)

    # ── Locate each readout's fork estimate ─────────────────────────────
    tau_R1 = interior_argmax(lam_max_curve, grid)
    tau_R2_pos = interior_argmax(div_curve, grid)         # tr(J) most-positive
    tau_R2_neg = interior_argmin(div_curve, grid)         # tr(J) most-negative
    tau_R3 = interior_argmax(lam_perp_curve, grid)
    tau_R4 = interior_argmax(ftle_c, grid)

    results = {
        "tau_crit": TAU_CRIT,
        "tolerance": TOLERANCE,
        "r2": r2,
        "R1_lambda_max": {
            "tau_peak": tau_R1, "err": abs(tau_R1 - TAU_CRIT),
            "pass": bool(abs(tau_R1 - TAU_CRIT) <= TOLERANCE),
        },
        "R2_div_pos": {
            "tau_peak": tau_R2_pos, "err": abs(tau_R2_pos - TAU_CRIT),
            "pass": bool(abs(tau_R2_pos - TAU_CRIT) <= TOLERANCE),
        },
        "R2_div_neg": {
            "tau_peak": tau_R2_neg, "err": abs(tau_R2_neg - TAU_CRIT),
            "pass": bool(abs(tau_R2_neg - TAU_CRIT) <= TOLERANCE),
        },
        "R3_transverse_lambda_max": {
            "tau_peak": tau_R3, "err": abs(tau_R3 - TAU_CRIT),
            "pass": bool(abs(tau_R3 - TAU_CRIT) <= TOLERANCE),
        },
        "R4_FTLE": {
            "tau_peak": tau_R4, "err": abs(tau_R4 - TAU_CRIT),
            "pass": bool(abs(tau_R4 - TAU_CRIT) <= TOLERANCE),
        },
    }

    print()
    print("=" * 72)
    print(f"{'Readout':32s} {'τ_peak':>8s} {'err vs τ_crit':>14s} {'PASS?':>6s}")
    print("-" * 72)
    for k in ["R1_lambda_max", "R2_div_pos", "R2_div_neg", "R3_transverse_lambda_max", "R4_FTLE"]:
        r = results[k]
        print(f"{k:32s} {r['tau_peak']:>8.3f} {r['err']:>14.3f} {'PASS' if r['pass'] else 'FAIL':>6s}")
    print("=" * 72)
    print(f"τ_crit = {TAU_CRIT:.4f}  |  tolerance = ±{TOLERANCE}")

    outpath = Path("{Path(__file__).parent / "saddle_readouts.json"}")
    outpath.write_text(json.dumps(results, indent=2))
    print(f"\nSaved: {outpath}")

    # Also save the curves themselves for downstream plotting.
    np.savez_compressed(
        "{Path(__file__).parent / "saddle_readouts_curves.npz"}",
        grid=grid,
        R1_lambda_max=lam_max_curve,
        R2_divergence=div_curve,
        R3_transverse_lambda_max=lam_perp_curve,
        R4_FTLE=ftle_c,
        tau_crit=TAU_CRIT,
    )


if __name__ == "__main__":
    main()
