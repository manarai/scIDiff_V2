"""
Temporal Jacobian tensor decomposition via semi-NMF.

Semi-NMF: J(t) ≈ Σ_k a_k(t) · A_k
  - activations a_k(t) ≥ 0  (interpretable: how active is archetype k at time t)
  - patterns    A_k  signed  (Jacobian operators, can have negative entries)

This matches the manuscript's "constrained factorization with non-negativity
on activations". SVD is kept as a fallback via jacobian_modes_svd().
"""
from __future__ import annotations

import numpy as np
import torch
from scipy.optimize import nnls


# ---------------------------------------------------------------------------
# Semi-NMF (primary)
# ---------------------------------------------------------------------------

def _nnls_rows(H: np.ndarray, M: np.ndarray) -> np.ndarray:
    """
    Solve W ≥ 0 minimising ||M - W @ H||_F row-wise via scipy NNLS.

    For each row m_t of M (shape D), solve min_{w_t ≥ 0} ||H.T @ w_t - m_t||_2,
    giving W of shape (T, K). scipy.optimize.nnls is the Lawson-Hanson active-set
    solver; each row is an independent QP, so the loop is embarrassingly parallel
    at the row level and cheap at typical (T, K) sizes.
    """
    T = M.shape[0]
    K = H.shape[0]
    A = H.T.astype(np.float64)  # (D, K)
    W = np.empty((T, K), dtype=np.float32)
    for t in range(T):
        w_t, _ = nnls(A, M[t].astype(np.float64))
        W[t] = w_t.astype(np.float32)
    return W


def _rescale(W: np.ndarray, H: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Resolve scale indeterminacy W → W D, H → D⁻¹ H by fixing ||H_k||_2 = 1
    for each archetype k. Preserves W @ H exactly (up to fp) and makes
    activation magnitudes comparable across runs and datasets.
    """
    h_norm = np.linalg.norm(H, axis=1)                # (K,)
    h_norm = np.where(h_norm < 1e-12, 1.0, h_norm)
    H_out = H / h_norm[:, None]
    W_out = W * h_norm[None, :]
    return W_out, H_out


def _seminmf(M: np.ndarray, rank: int, max_iter: int = 500,
             tol: float = 1e-4, seed: int = 0,
             tv_lambda: float = 0.0) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Semi-NMF via alternating least squares.
    M (T, D) ≈ W (T, K) @ H (K, D)  where W ≥ 0, H unconstrained.

    W-update uses column-wise NNLS (scipy.optimize.nnls, Lawson-Hanson) so the
    Frobenius objective is monotone non-increasing across iterations. The prior
    implementation clipped the unconstrained ridge solution, which does not
    solve the constrained subproblem and can increase the objective.

    Parameters
    ----------
    tv_lambda : float, optional (default 0.0)
        Total-variation smoothness penalty on the temporal derivative of W.
        When > 0, the W update switches from closed-form NNLS to a projected
        proximal-gradient inner loop that penalises Σ_t ||W[t+1] - W[t]||_1
        column-wise. This gives the algorithm a preference for smooth (and,
        at strong penalty, coherent) activation profiles — useful for
        recovering concurrent-Gaussian activations that semi-NMF otherwise
        splits into a spread pattern because the unregularised objective is
        indifferent between concurrent and spread solutions. See
        :func:`jacobian_modes` for guidance on when to enable this.
    """
    rng = np.random.default_rng(seed)
    T, D = M.shape

    W = np.abs(rng.standard_normal((T, rank)).astype(np.float32)) + 1e-4
    W /= W.sum(axis=0, keepdims=True) + 1e-8

    prev_err = np.inf
    err = np.inf
    for _ in range(max_iter):
        # H update — unconstrained least squares.
        WtW = W.T @ W + 1e-6 * np.eye(rank)
        H = np.linalg.solve(WtW, W.T @ M)  # (K, D)

        # W update.
        if tv_lambda <= 0.0:
            W = _nnls_rows(H, M)
        else:
            W = _tv_smoothed_W_update(W, H, M, tv_lambda)

        err = float(np.linalg.norm(M - W @ H, "fro"))
        if abs(prev_err - err) / (prev_err + 1e-8) < tol:
            break
        prev_err = err

    W, H = _rescale(W, H)
    return W, H, err


def _tv_smoothed_W_update(W: np.ndarray, H: np.ndarray, M: np.ndarray,
                          tv_lambda: float, n_inner: int = 40) -> np.ndarray:
    """
    ISTA-style W update for the objective
        ||M - W H||_F^2 + tv_lambda * Σ_k Σ_t |W[t+1, k] - W[t, k]|,
    with W ≥ 0.

    Gradient of the smooth part w.r.t. W is 2 (W H - M) H^T. The
    proximal operator of the TV term separates per column and reduces to
    soft-thresholding of the first differences (cumulative sum after
    thresholding). Non-negativity is enforced by projection after the prox.
    """
    rank = W.shape[1]
    norm_H = float(np.linalg.norm(H, 2))
    step = 1.0 / (2.0 * norm_H ** 2 + 1e-8)
    thresh = tv_lambda * step
    W_out = W.astype(np.float32).copy()
    for _ in range(n_inner):
        grad = 2.0 * (W_out @ H - M) @ H.T
        W_prox = W_out - step * grad
        for k in range(rank):
            col = W_prox[:, k]
            diffs = np.diff(col)
            sign = np.sign(diffs)
            mag = np.maximum(np.abs(diffs) - thresh, 0.0)
            new_diffs = sign * mag
            col_new = np.empty_like(col)
            col_new[0] = col[0]
            col_new[1:] = col[0] + np.cumsum(new_diffs)
            W_prox[:, k] = col_new
        W_out = np.maximum(W_prox, 0.0).astype(np.float32)
    return W_out


def jacobian_modes(
    J_tensor: torch.Tensor,
    rank: int = 5,
    n_restarts: int = 5,
    seed: int = 0,
    tv_lambda: float = 0.0,
) -> tuple[torch.Tensor, torch.Tensor, float]:
    """
    Decompose a temporal Jacobian tensor into recurrent operator archetypes.

    Parameters
    ----------
    J_tensor   : (T, D, D) — stacked Jacobians across pseudotime windows.
    rank       : Number of archetypes K.
    n_restarts : ALS restarts (best reconstruction kept).
    seed       : Base random seed.
    tv_lambda  : Total-variation smoothness penalty on temporal derivatives of
        the activation profiles (rows-of-time in the W factor). Default 0.0
        preserves the unregularised behaviour used for the main analyses.
        Enable (typical range 0.01–1.0) when the biology of interest involves
        concurrent (co-timed) archetype activations: the unregularised semi-NMF
        objective is indifferent between a "concurrent peaks" solution and a
        "spread peaks" solution that reconstructs the same tensor equally well,
        and picks the spread one, so concurrent coordination gets under-called
        without this prior. Sequential-handoff structure does NOT require this
        penalty; leave at 0.0 for those analyses. On the Figure 2 System 2
        synthetic benchmark, enabling tv_lambda moves matched
        concurrent-peak-pair prevalence from 0.05 toward the ground-truth 1.0
        (at the cost of higher reconstruction error at strong penalties).

    Returns
    -------
    patterns    : (K, D, D) — archetype Jacobian matrices (signed).
    activations : (T, K)    — non-negative temporal activation profiles.
    error       : float     — final Frobenius reconstruction error.
    """
    if J_tensor.dim() == 4:
        J_tensor = J_tensor.mean(dim=1)

    T, d1, d2 = J_tensor.shape
    M = J_tensor.reshape(T, d1 * d2).numpy().astype(np.float32)

    best_W, best_H, best_err = None, None, np.inf
    for s in range(n_restarts):
        W, H, err = _seminmf(M, rank=rank, seed=seed + s, tv_lambda=tv_lambda)
        if err < best_err:
            best_W, best_H, best_err = W, H, err

    patterns    = torch.tensor(np.asarray(best_H), dtype=torch.float32).reshape(rank, d1, d2)
    activations = torch.tensor(np.asarray(best_W), dtype=torch.float32)
    return patterns, activations, best_err


# ---------------------------------------------------------------------------
# SVD fallback (kept for backward compatibility)
# ---------------------------------------------------------------------------

def jacobian_modes_svd(
    J_tensor: torch.Tensor,
    rank: int = 5,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """SVD-based decomposition (activations can be negative). Legacy API."""
    if J_tensor.dim() == 4:
        J_tensor = J_tensor.mean(dim=1)
    T, d1, d2 = J_tensor.shape
    M          = J_tensor.reshape(T, d1 * d2)
    U, S, Vh   = torch.linalg.svd(M, full_matrices=False)
    patterns    = Vh[:rank].reshape(rank, d1, d2)
    activations = U[:, :rank]
    return patterns, activations, S[:rank]
