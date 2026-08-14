"""
Tests for the concurrent-activation recovery gap in Figure 2 System 2.

Hypothesis (from the critique): the 0.06 matched vs 1.00 GT concurrent-peak
prevalence is NOT a data-identifiability limit — semi-NMF has no temporal
model, so concurrent-Gaussian activations that are near-collinear in the
unfolded representation get split into a spread pattern by the unconstrained
factorization. If the algorithm has no preference between "all peaks near
t=0.5" and "peaks spread evenly", it will pick spread.

Three cheap tests, ordered by cost:

  T1  Reconstruction-residual comparison. Compute the Frobenius residual at
      the ground-truth (W_true, H_true) vs the semi-NMF-recovered (W, H).
      If they are comparable, the objective cannot distinguish concurrent
      from spread, and the fix is a prior; if the recovered has lower error,
      the recovered activations are genuinely a better fit by MSE and the
      "recovery gap" is a metric-vs-objective mismatch, not a decomposition
      failure.

  T2  K sweep. If cross-validated K selects K > true K, extra components
      split concurrent programs mechanically. Sweep K ∈ {2, 3, 4, 5, 6, 7}
      on the true-K=5 system and measure concurrent_prevalence at each K.

  T3  SVD fallback. jacobian_modes_svd already exists and has an implicit
      orthogonality prior on left factors. Compare on the same tensor.

If T1 + T2 + T3 confirm the diagnosis, T4 is the fix:
  T4  TV-smoothed semi-NMF: add a fused-lasso penalty λ Σ|H[t+1]-H[t]| on
      rows of H so the algorithm has a preference for smooth (and, at high
      enough λ, coherent) activation profiles.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np
import torch
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import linear_sum_assignment, nnls
from scipy.stats import kendalltau

REPO = Path("/Users/terooatt/Downloads/scJDO")
sys.path.insert(0, str(REPO))
from scjdo.archetypes.decompose import _seminmf, jacobian_modes, jacobian_modes_svd

# ── Figure 2 System 2 config (from Figures_notebook/Figure2_synthetic_peak_timing_null_analysis.ipynb) ──
METHODS_N_WINDOWS = 100
METHODS_N_PCS = 50
RANK_TRUE = 5
N_SEEDS = 8
GLOBAL_SEED = 20260515
TENSOR_NOISE = 0.010
CONCURRENT_MAX_PEAK_GAP = 0.045


def stable_seed(*items, base=GLOBAL_SEED):
    payload = "|".join(map(str, items)).encode("utf-8")
    digest = hashlib.sha256(payload).hexdigest()[:12]
    return (base + int(digest, 16)) % (2**32 - 1)


def normalize01(x, axis=0):
    x = np.asarray(x, dtype=np.float32)
    return (x - x.min(axis=axis, keepdims=True)) / (
        x.max(axis=axis, keepdims=True) - x.min(axis=axis, keepdims=True) + 1e-8
    )


def gaussian_bump(t, center, width):
    return np.exp(-0.5 * ((t - center) / width) ** 2)


def system2_profiles(seed_index):
    """System 2: concurrent Gaussians with varying overlap widths."""
    rng = np.random.default_rng(stable_seed("system_2_concurrent_varying_overlap", seed_index, "profiles"))
    t = np.linspace(0, 1, METHODS_N_WINDOWS, dtype=np.float32)
    center = 0.52 + rng.normal(0, 0.006)
    widths = np.array([0.075, 0.105, 0.140, 0.180, 0.220]) * rng.uniform(0.88, 1.12, RANK_TRUE)
    A = []
    for width in widths:
        local_center = center + rng.normal(0, 0.008)
        prof = gaussian_bump(t, local_center, width)
        prof += 0.020 * gaussian_bump(t, local_center + rng.normal(0, 0.006), width * 0.60)
        prof += 0.006 * gaussian_filter1d(rng.standard_normal(METHODS_N_WINDOWS), sigma=1.5)
        A.append(prof)
    return t, normalize01(np.stack(A, axis=1), axis=0).astype(np.float32)


def make_operator_basis(seed, d=METHODS_N_PCS, rank=RANK_TRUE):
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.standard_normal((d, d)))
    basis = []
    for k in range(rank):
        M = np.zeros((d, d), dtype=np.float32)
        lo = 2 * k
        M[lo:lo+2, lo:lo+2] = np.array(
            [[(-1)**k * (0.7 + 0.1*k), -0.6], [0.6, -0.4 - 0.05*k]], dtype=np.float32
        )
        if lo + 4 <= d:
            M[lo+2:lo+4, lo+2:lo+4] = np.diag([0.5 + 0.05*k, -0.3 - 0.05*k])
        dense = Q @ M @ Q.T
        basis.append((dense / (np.linalg.norm(dense, "fro") + 1e-8)).astype(np.float32))
    return np.stack(basis).astype(np.float32)


def make_tensor(A, seed):
    rng = np.random.default_rng(seed)
    B = make_operator_basis(seed=seed + 17)
    J = np.einsum("tk,kij->tij", A, B).astype(np.float32)
    eps = rng.standard_normal(J.shape).astype(np.float32)
    eps = gaussian_filter1d(eps, sigma=1.0, axis=0)
    J = J + TENSOR_NOISE * eps / (eps.std() + 1e-8) * (J.std() + 1e-8)
    return J.astype(np.float32), B


def concurrent_prevalence(A, t, gap=CONCURRENT_MAX_PEAK_GAP):
    """Fraction of archetype pairs whose peaks are within `gap` of each other."""
    peaks = t[np.argmax(A, axis=0)]
    K = A.shape[1]
    if K < 2:
        return 0.0
    pairs = 0
    hit = 0
    for i in range(K):
        for j in range(i + 1, K):
            pairs += 1
            if abs(peaks[i] - peaks[j]) <= gap:
                hit += 1
    return hit / pairs


def match_columns(A_gt, A_rec):
    """Hungarian match on activation cosine (used only for T1's residual comparison)."""
    Kg = A_gt.shape[1]
    Kr = A_rec.shape[1]
    K = min(Kg, Kr)
    ag = A_gt[:, :K] / (np.linalg.norm(A_gt[:, :K], axis=0, keepdims=True) + 1e-8)
    ar = A_rec[:, :K] / (np.linalg.norm(A_rec[:, :K], axis=0, keepdims=True) + 1e-8)
    sim = ag.T @ ar
    row, col = linear_sum_assignment(-sim)
    perm = np.full(K, -1, dtype=int)
    for r, c in zip(row, col):
        perm[r] = c
    return perm


# ── T4: TV-smoothed semi-NMF ────────────────────────────────────────────
# The unfolded matrix M has shape (T, D*D). Semi-NMF factorises M ≈ W H,
# W ≥ 0, H signed. We already know the row structure of M is time-ordered;
# adding a total-variation penalty λ Σ_{t} ||W[t+1] - W[t]||_1 on the W
# activations (which are the time-varying profiles) encodes that.

def _tv_seminmf(M, rank, tv_lambda, max_iter=300, tol=1e-4, seed=0):
    """Semi-NMF with a total-variation penalty on the temporal axis of W.

    W update: still column-wise NNLS on the augmented system that stacks the
    reconstruction and the TV difference operator D_t. Concretely:

        min_{W >= 0} || [ H^T ; sqrt(lambda) * D_t ] W_t - [ M_t ; 0 ] ||^2

    is not separable across t (D_t couples adjacent rows). Two practical
    options: (a) proximal gradient descent on the W step, (b) block-coordinate
    with augmented linear system. We use (a) — projected gradient + soft
    threshold — which is short to write and easy to verify converges on this
    small (T=100, K<=7) size.

    H update remains the closed-form LS solution (H unconstrained).
    """
    rng = np.random.default_rng(seed)
    T, D = M.shape
    W = np.abs(rng.standard_normal((T, rank)).astype(np.float32)) + 1e-4
    W /= W.sum(axis=0, keepdims=True) + 1e-8

    # First-difference operator D_t of shape (T-1, T) so D_t @ W has shape (T-1, K).
    prev_err = np.inf
    err = np.inf
    for it in range(max_iter):
        # H update — closed form.
        WtW = W.T @ W + 1e-6 * np.eye(rank)
        H = np.linalg.solve(WtW, W.T @ M)

        # W update — proximal gradient descent with non-negativity projection
        # on the objective ||M - W H||_F^2 + lambda * ||D_t W||_1.
        # Grad of Frobenius part w.r.t. W: 2 (W H - M) H^T
        # ISTA step with step size = 1 / (2 * ||H||_2^2)
        norm_H = float(np.linalg.norm(H, 2))
        step = 1.0 / (2.0 * norm_H ** 2 + 1e-8)
        for _inner in range(40):
            grad = 2.0 * (W @ H - M) @ H.T
            W_prox = W - step * grad
            # Soft-threshold the temporal differences on each column.
            for k in range(rank):
                col = W_prox[:, k]
                diffs = np.diff(col)
                sign = np.sign(diffs)
                mag = np.maximum(np.abs(diffs) - tv_lambda * step, 0.0)
                new_diffs = sign * mag
                col_new = np.empty_like(col)
                col_new[0] = col[0]
                col_new[1:] = col[0] + np.cumsum(new_diffs)
                W_prox[:, k] = col_new
            W_prox = np.maximum(W_prox, 0.0)
            W = W_prox

        err = float(np.linalg.norm(M - W @ H, "fro"))
        if abs(prev_err - err) / (prev_err + 1e-8) < tol:
            break
        prev_err = err

    # Rescale so rows of H have unit norm.
    h_norm = np.linalg.norm(H, axis=1)
    h_norm = np.where(h_norm < 1e-12, 1.0, h_norm)
    H = H / h_norm[:, None]
    W = W * h_norm[None, :]
    return W.astype(np.float32), H.astype(np.float32), err


# ── Runners ─────────────────────────────────────────────────────────────

def run_one_seed(seed_index):
    t, A_gt = system2_profiles(seed_index)
    J, B = make_tensor(A_gt, seed=stable_seed("system_2_concurrent_varying_overlap", seed_index, "tensor"))
    return t, A_gt, J


def T1_residual_at_gt_vs_recovered(t, A_gt, J):
    """Reconstruct from GT (W=A_gt, H=B.reshape) vs NMF-recovered; compare residuals."""
    T = A_gt.shape[0]
    d = J.shape[1]
    M = J.reshape(T, d * d).astype(np.float32)
    # GT reconstruction: J_true = A_gt @ B_reshaped
    B_shape = M @ np.linalg.pinv(A_gt @ A_gt.T @ A_gt.T) if False else None  # not needed here
    # The true B is what we passed into make_tensor. Instead of trying to recover it,
    # just fit H from A_gt: H = pinv(A_gt) M, which is the LS-best H given W = A_gt.
    W_gt = A_gt
    H_from_gt_W = np.linalg.pinv(W_gt) @ M
    res_gt = float(np.linalg.norm(M - W_gt @ H_from_gt_W, "fro"))

    # NMF-recovered
    patterns_nmf, W_nmf, err_nmf = jacobian_modes(
        torch.from_numpy(J), rank=RANK_TRUE, n_restarts=5, seed=0
    )
    W_nmf_np = W_nmf.numpy()
    H_nmf_np = patterns_nmf.reshape(RANK_TRUE, -1).numpy()
    res_nmf = float(np.linalg.norm(M - W_nmf_np @ H_nmf_np, "fro"))
    return res_gt, res_nmf, W_nmf_np


def T2_K_sweep(J, t, A_gt, K_values=(2, 3, 4, 5, 6, 7)):
    out = {}
    for K in K_values:
        patterns, activations, err = jacobian_modes(
            torch.from_numpy(J), rank=K, n_restarts=5, seed=0
        )
        A_rec = normalize01(activations.numpy(), axis=0)
        # Match columns (only used for coherent labelling — concurrent metric is permutation-invariant)
        cp = concurrent_prevalence(A_rec, t)
        peaks = t[np.argmax(A_rec, axis=0)]
        out[K] = {
            "concurrent_prevalence": float(cp),
            "peak_spread": float(peaks.max() - peaks.min()),
            "reconstruction_err": float(err),
        }
    return out


def T3_svd_vs_seminmf(J, t):
    # SVD
    patterns_svd, activations_svd, S = jacobian_modes_svd(
        torch.from_numpy(J), rank=RANK_TRUE
    )
    A_svd = normalize01(activations_svd.numpy(), axis=0)
    cp_svd = concurrent_prevalence(A_svd, t)

    # semi-NMF for comparison
    patterns_nmf, activations_nmf, err_nmf = jacobian_modes(
        torch.from_numpy(J), rank=RANK_TRUE, n_restarts=5, seed=0
    )
    A_nmf = normalize01(activations_nmf.numpy(), axis=0)
    cp_nmf = concurrent_prevalence(A_nmf, t)

    return {"svd": float(cp_svd), "seminmf": float(cp_nmf)}


def T4_tv_seminmf(J, t, tv_values=(0.0, 0.01, 0.1, 1.0, 10.0)):
    T = J.shape[0]
    d = J.shape[1]
    M = J.reshape(T, d * d).astype(np.float32)
    out = {}
    for tv in tv_values:
        W, H, err = _tv_seminmf(M, rank=RANK_TRUE, tv_lambda=tv, seed=0)
        A = normalize01(W, axis=0)
        cp = concurrent_prevalence(A, t)
        out[tv] = {
            "concurrent_prevalence": float(cp),
            "reconstruction_err": float(err),
        }
    return out


def main():
    print("=" * 76)
    print("Concurrent-activation recovery diagnostics — Figure 2 System 2")
    print("=" * 76)

    # Aggregate over seeds
    t1_res_gt, t1_res_nmf = [], []
    t2_by_seed = []
    t3_by_seed = []
    t4_by_seed = []
    gt_cp = []

    for si in range(N_SEEDS):
        t, A_gt, J = run_one_seed(si)
        gt_cp.append(concurrent_prevalence(A_gt, t))

        r_gt, r_nmf, _ = T1_residual_at_gt_vs_recovered(t, A_gt, J)
        t1_res_gt.append(r_gt)
        t1_res_nmf.append(r_nmf)

        t2_by_seed.append(T2_K_sweep(J, t, A_gt))
        t3_by_seed.append(T3_svd_vs_seminmf(J, t))
        t4_by_seed.append(T4_tv_seminmf(J, t))

    print(f"\nGT concurrent_prevalence (mean over {N_SEEDS} seeds): {np.mean(gt_cp):.3f}")

    # T1
    print("\n--- T1  Residual at ground-truth W vs recovered W ---")
    print(f"  Frobenius residual at W = A_gt (H = LS-best given W): "
          f"{np.mean(t1_res_gt):.4f} ± {np.std(t1_res_gt):.4f}")
    print(f"  Frobenius residual at W = semi-NMF recovered:         "
          f"{np.mean(t1_res_nmf):.4f} ± {np.std(t1_res_nmf):.4f}")
    ratio = np.mean(t1_res_nmf) / np.mean(t1_res_gt)
    print(f"  Ratio (recovered / GT): {ratio:.3f}")
    if ratio < 1.0:
        print(f"  → Recovered residual is LOWER — the objective genuinely prefers the "
              f"spread solution over the true one. Adding a temporal prior is the "
              f"principled fix (T4).")
    elif ratio > 1.05:
        print(f"  → Recovered residual is HIGHER — the algorithm is finding a suboptimum, "
              f"not a fundamentally-different objective preference. More restarts or a "
              f"better initialiser may help.")
    else:
        print(f"  → Residuals are comparable — the objective is nearly indifferent, and a "
              f"prior is the appropriate lever.")

    # T2
    print("\n--- T2  K sweep on true-K=5 system ---")
    print(f"  {'K':>4}  {'concurrent_prev':>16}  {'peak_spread':>12}  {'recon_err':>10}")
    for K in (2, 3, 4, 5, 6, 7):
        vals_cp = [seed[K]["concurrent_prevalence"] for seed in t2_by_seed]
        vals_ps = [seed[K]["peak_spread"] for seed in t2_by_seed]
        vals_re = [seed[K]["reconstruction_err"] for seed in t2_by_seed]
        print(f"  {K:>4}  {np.mean(vals_cp):>16.3f}  {np.mean(vals_ps):>12.3f}  {np.mean(vals_re):>10.4f}")

    # T3
    print("\n--- T3  SVD vs semi-NMF (K = 5) ---")
    for method in ("svd", "seminmf"):
        vals = [seed[method] for seed in t3_by_seed]
        print(f"  {method:>10}: concurrent_prevalence = {np.mean(vals):.3f} ± {np.std(vals):.3f}")

    # T4
    print("\n--- T4  TV-smoothed semi-NMF (K = 5, sweep TV weight) ---")
    print(f"  {'lambda':>8}  {'concurrent_prev':>16}  {'recon_err':>10}")
    for tv in (0.0, 0.01, 0.1, 1.0, 10.0):
        vals_cp = [seed[tv]["concurrent_prevalence"] for seed in t4_by_seed]
        vals_re = [seed[tv]["reconstruction_err"] for seed in t4_by_seed]
        print(f"  {tv:>8}  {np.mean(vals_cp):>16.3f}  {np.mean(vals_re):>10.4f}")

    print("\n" + "=" * 76)
    print("Verdict:")
    print("  Baseline concurrent_prevalence (K=5, semi-NMF, no TV): "
          f"{np.mean([seed[5]['concurrent_prevalence'] for seed in t2_by_seed]):.3f} (GT = 1.00)")
    print("  Best TV concurrent_prevalence:                          "
          f"{max(np.mean([seed[tv]['concurrent_prevalence'] for seed in t4_by_seed]) for tv in (0.0, 0.01, 0.1, 1.0, 10.0)):.3f}")
    print("=" * 76)


if __name__ == "__main__":
    main()
