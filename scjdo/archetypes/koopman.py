"""
Windowed, nonautonomous Koopman decomposition of the temporal Jacobian
tensor produced by scJDO.

Rationale
---------
scJDO learns a drift field f(x, t) and its Jacobian J(t) = ∂f/∂x along
pseudotime. The temporal tensor J(t) ∈ R^{T x D x D} captures how local
sensitivity evolves. Semi-NMF factorises this tensor into non-negative
activations of signed operator archetypes; Koopman gives a complementary
spectral lens by treating the vectorised Jacobian sequence as a
trajectory of observables and fitting a linear operator on it. The
eigenvalues of that operator are growth/decay rates and oscillation
frequencies; the eigenvectors, un-flattened, act as operator "modes".

Model
-----
Let y(t) = vec(J(t)) ∈ R^{D^2} and z(t) = P_r y(t) ∈ R^r its rank-r
projection (P_r from SVD of Y = [y(t_1); ...; y(t_T)]). EDMD with
identity observables (i.e. DMD in the reduced space) solves, for a
sliding window W_τ of consecutive snapshots,

    K_τ = argmin_K  Σ_{i ∈ W_τ}  || z(t_{i+1}) - K z(t_i) ||^2,

via least squares. Eigendecompose K_τ v_k = λ_k v_k. With average
window pseudotime step Δτ, the continuous-time rate of mode k is
μ_k = log(λ_k) / Δτ, so Re(μ_k) > 0 marks growth (instability) and
|Im(μ_k)| > 0 marks oscillation with angular frequency |Im(μ_k)|.

Windowing is per-grid-point because single-cell trajectories are
typically nonstationary; a single global K blurs branch- and
regime-specific behaviour. A ``mode='global'`` fallback is offered as a
baseline.

Output shape parity with semi-NMF
---------------------------------
:func:`koopman_modes` returns (patterns, activations, error) with the
same shapes as :func:`scjdo.archetypes.decompose.jacobian_modes`, so
downstream code (``_gene_scores``, per-window instability projection,
etc.) does not need to change. Complex conjugate Koopman
eigenvectors are split into real / imaginary parts, which keeps the
returned pattern matrices real-valued.

Caveats specific to snapshot scRNA-seq
--------------------------------------
- Pseudotime is an inferred surrogate for physical time; the operator
  here summarises the inferred progression, not the intrinsic flow.
- Snapshot data cannot unambiguously identify oscillatory vs.
  non-oscillatory flows — treat |Im(μ_k)| > 0 modes as candidate
  oscillations, not as physically resolved cycles.
- Boundary windows suffer spectral pollution; interior grid points are
  more trustworthy.
"""
from __future__ import annotations

import warnings
from typing import Optional

import numpy as np
import torch


# ---------------------------------------------------------------------------
# Low-level EDMD (identity observables — i.e. reduced DMD)
# ---------------------------------------------------------------------------

def _fit_edmd(Z: np.ndarray, Zp: np.ndarray, ridge: float = 0.0) -> np.ndarray:
    """
    Fit K minimising || Zp - Z K^T ||_F on snapshot pairs.

    Z, Zp have shape (n_pairs, r): each row is a reduced snapshot. The
    least-squares operator is K = (Z^+ Zp)^T. A small Tikhonov ridge
    stabilises the solve when the window is short or nearly rank
    deficient.
    """
    if ridge > 0.0:
        A = Z.T @ Z + ridge * np.eye(Z.shape[1], dtype=Z.dtype)
        B = Z.T @ Zp
        K_T = np.linalg.solve(A, B)          # r x r
    else:
        K_T, *_ = np.linalg.lstsq(Z, Zp, rcond=None)
    return K_T.T.astype(np.float64)          # K: r x r


def _reduce_trajectory(Y: np.ndarray, rank: int,
                       whiten: bool = False) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    SVD reduction of the trajectory matrix Y (T x D²) to rank r.

    Two metric choices; the diagnostics geometry (Henrici, reactivity,
    transient gain, eigvec conditioning) depends on which is used.
    Eigenvalues of the fitted operator do not — the two Ks are similar,
    so their spectra coincide.

    Parameters
    ----------
    whiten : bool
        - ``False`` (default) — ship behaviour. ``Z = Y V_r = U_r Σ_r``,
          carrying the singular-value scaling. Directions with larger
          singular values dominate the least-squares fit; the reduced
          space has a non-neutral metric.
        - ``True`` — orthonormalise the trajectory in the reduced basis.
          ``Z = U_r`` with orthonormal columns; the fit is balanced
          across modes. Under exact least squares the two operators are
          similar (same eigenvalues, different eigenvectors); with
          ``ridge > 0`` the regularisation acts differently in the two
          metrics, so eigenvalues shift slightly as well.

    Returns
    -------
    Z    : (T, r)    reduced trajectory (rows are per-time coordinates).
    V_r  : (D², r)   right singular vectors — always the same regardless
                     of ``whiten``; used as an orthonormal basis for the
                     vectorised operator.
    S_r  : (r,)      singular values.
    lift : (D², r)   the matrix that lifts an r-vector expressed in the
                     *current* reduced coordinates back to D². Equal to
                     ``V_r`` when ``whiten=False`` and ``V_r Σ_r`` when
                     ``whiten=True``.
    """
    U, S, Vt = np.linalg.svd(Y, full_matrices=False)
    r = min(rank, U.shape[1])
    V_r = Vt[:r].T                            # (D², r)
    S_r = S[:r]
    if whiten:
        Z = U[:, :r].astype(np.float64)                       # orthonormal columns
        lift = (V_r * S_r[None, :]).astype(np.float64)        # V_r Σ_r
    else:
        Z = (Y @ V_r).astype(np.float64)                      # U_r Σ_r
        lift = V_r.astype(np.float64)                         # V_r
    return Z, V_r.astype(np.float64), S_r.astype(np.float64), lift


def _windows(T: int, half: int) -> list[np.ndarray]:
    """Sliding [t-half, t+half] index arrays with edge clipping."""
    return [np.arange(max(0, t - half), min(T, t + half + 1)) for t in range(T)]


# ---------------------------------------------------------------------------
# Discrete-time operator geometry — branch-free descriptors
# ---------------------------------------------------------------------------

def _matrix_geometry(K: np.ndarray, n_max: int = 30) -> dict:
    """
    Geometry descriptors computed directly on the discrete-time operator
    K, without taking a logarithm — so they are immune to the
    ``log λ = log λ + 2π i k`` branch ambiguity that constrains the
    frequency descriptors. In the snapshot-scRNA-seq setting where the
    pseudotime step Δτ is coarse relative to the fastest true dynamics,
    these tend to be the better-posed half of the spectral output.

    Descriptors
    -----------
    henrici        : ν = sqrt(||K||_F² − Σ|λ_i|²) / ||K||_F, in [0, 1].
                     Departure-from-normality index — 0 iff K is normal.
    reactivity     : σ_max(K) / ρ(K), where σ_max is the operator norm
                     and ρ the spectral radius. Always ≥ 1; equality iff
                     K is normal. Reads as "how much can a state be
                     amplified in one step relative to its
                     asymptotic rate".
    transient_gain : max_{n ∈ [1, n_max]} ||K^n||_2. Captures short-time
                     amplification that asymptotic eigenvalues can miss
                     when K is highly non-normal.
    eigvec_cond    : condition number of the (complex) eigenvector
                     matrix V. Blows up when K is nearly defective; a
                     large value warns that the eigen-based
                     interpretation is fragile.

    All four depend on the reduced-space metric (i.e. on the ``whiten``
    choice at :func:`_reduce_trajectory`). Document which was used.
    """
    K = np.asarray(K, dtype=np.float64)
    r = K.shape[0]
    frob2 = float(np.sum(K * K))
    frob = float(np.sqrt(frob2))
    eigs = np.linalg.eigvals(K)
    sum_lam2 = float(np.sum(np.abs(eigs) ** 2))
    henrici_abs = float(np.sqrt(max(0.0, frob2 - sum_lam2)))
    henrici = henrici_abs / (frob + 1e-12)

    sigma_max = float(np.linalg.norm(K, ord=2))
    rho = float(np.max(np.abs(eigs))) if eigs.size else 0.0
    reactivity = sigma_max / (rho + 1e-12)

    P = np.eye(r, dtype=np.float64)
    gain = 0.0
    for _ in range(int(n_max)):
        P = K @ P
        gn = float(np.linalg.norm(P, ord=2))
        if gn > gain:
            gain = gn
        # If the operator is contracting and has dropped far below
        # numerical noise, further powers only add roundoff.
        if rho < 1.0 and gn < 1e-12:
            break

    try:
        _, V = np.linalg.eig(K)
        kappa = float(np.linalg.cond(V))
    except np.linalg.LinAlgError:
        kappa = float("inf")

    return dict(
        henrici=henrici,
        reactivity=reactivity,
        transient_gain=gain,
        eigvec_cond=kappa,
    )


# ---------------------------------------------------------------------------
# Complex-mode → real-mode splitting
# ---------------------------------------------------------------------------

def _split_complex_pairs(eigvals: np.ndarray, eigvecs: np.ndarray,
                         tol: float = 1e-8) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Sort eigenpairs by |λ| descending and split complex conjugates into
    (real part, imaginary part) so the resulting basis is real-valued.

    Returns (lam_all, V_real, is_imag_component) where V_real columns are
    real vectors (Re v for the real part of a complex pair, Im v for its
    imaginary partner, or v itself for a real eigenvalue), and
    is_imag_component[i] tells whether column i was the imaginary partner
    (useful for downstream sign / phase handling).
    """
    order = np.argsort(-np.abs(eigvals))
    lam_sorted = eigvals[order]
    V_sorted = eigvecs[:, order]
    r = V_sorted.shape[0]

    lam_out: list[complex] = []
    V_out: list[np.ndarray] = []
    is_imag: list[bool] = []
    taken = np.zeros(len(lam_sorted), dtype=bool)
    for i in range(len(lam_sorted)):
        if taken[i]:
            continue
        taken[i] = True
        lam_i = lam_sorted[i]
        v_i = V_sorted[:, i]
        if abs(lam_i.imag) < tol:
            lam_out.append(complex(lam_i.real, 0.0))
            V_out.append(v_i.real.astype(np.float64))
            is_imag.append(False)
            continue
        # Find the conjugate partner.
        partner = -1
        for j in range(i + 1, len(lam_sorted)):
            if taken[j]:
                continue
            if abs(lam_sorted[j] - np.conj(lam_i)) < 1e-6 * max(1.0, abs(lam_i)):
                partner = j
                break
        if partner < 0:
            # No partner (rare, numerical) — fall back to the real projection.
            lam_out.append(complex(lam_i.real, 0.0))
            V_out.append(v_i.real.astype(np.float64))
            is_imag.append(False)
        else:
            taken[partner] = True
            lam_out.extend([lam_i, np.conj(lam_i)])
            V_out.extend([v_i.real.astype(np.float64),
                          v_i.imag.astype(np.float64)])
            is_imag.extend([False, True])

    V_real = np.stack(V_out, axis=1) if V_out else np.zeros((r, 0))
    return np.asarray(lam_out), V_real.astype(np.float64), np.asarray(is_imag)


# ---------------------------------------------------------------------------
# Local (windowed / nonautonomous) Koopman
# ---------------------------------------------------------------------------

def _local_koopman(
    Z: np.ndarray,
    grid: np.ndarray,
    window_half: int,
    ridge: float,
) -> tuple[list[np.ndarray], np.ndarray]:
    """
    Fit one local operator per grid point via sliding-window EDMD.

    Returns (K_seq, dt_seq) with K_seq[t] the r x r local operator and
    dt_seq[t] the mean pseudotime spacing in that window (used to
    convert discrete eigenvalues into continuous-time rates).
    """
    T = Z.shape[0]
    Zpair_lhs = Z[:-1]
    Zpair_rhs = Z[1:]
    win = _windows(T - 1, window_half)
    K_seq: list[np.ndarray] = []
    dt_seq = np.zeros(T, dtype=np.float64)
    for t in range(T):
        # Use windows on the *pair index* (there are T-1 pairs).
        wt = win[min(t, T - 2)]
        K_t = _fit_edmd(Zpair_lhs[wt], Zpair_rhs[wt], ridge=ridge)
        K_seq.append(K_t)
        span = grid[wt + 1] - grid[wt]
        dt_seq[t] = float(np.mean(span)) if len(span) else float(np.median(np.diff(grid)))
    return K_seq, dt_seq


def _global_koopman(Z: np.ndarray, ridge: float) -> np.ndarray:
    """Single global operator baseline (autonomous, stationary K)."""
    return _fit_edmd(Z[:-1], Z[1:], ridge=ridge)


# ---------------------------------------------------------------------------
# Public API — matches jacobian_modes signature
# ---------------------------------------------------------------------------

def koopman_modes(
    J_tensor: torch.Tensor,
    rank: int = 5,
    *,
    reduced_rank: Optional[int] = None,
    mode: str = "local",
    window_half: int = 8,
    ridge: float = 1e-4,
    grid: Optional[np.ndarray] = None,
    smooth_sigma: float = 1.0,
    whiten: bool = False,
    geometry_n_max: int = 30,
    return_diagnostics: bool = False,
) -> tuple[torch.Tensor, torch.Tensor, float] | tuple[torch.Tensor, torch.Tensor, float, dict]:
    """
    Windowed Koopman decomposition of a temporal Jacobian tensor.

    Parameters
    ----------
    J_tensor
        (T, D, D) stacked Jacobians. (T, N, D, D) is also accepted and is
        averaged over the second axis first, matching jacobian_modes.
    rank
        Number of Koopman modes ``K`` returned in the (K, D, D) patterns
        tensor. Complex conjugate pairs count as two modes (real + imag).
    reduced_rank
        Rank ``r`` of the SVD reduction on the vectorised trajectory
        before EDMD. Defaults to ``min(T-1, 4*rank+8)``. Keep r small
        relative to T so the local least-squares problem is well posed.
    mode
        ``'local'`` (default) fits one operator per grid point on a
        sliding window (nonautonomous Koopman). ``'global'`` fits a
        single stationary operator on all snapshot pairs (baseline).
    window_half
        Half-width of the sliding window in grid steps. Effective window
        length is 2*window_half+1. Larger = smoother spectrum, worse
        temporal resolution.
    ridge
        Tikhonov ridge on the reduced least-squares solve. Small (1e-6
        to 1e-3) is usually enough; increase when windows are short.
    grid
        Optional (T,) pseudotime grid. Used only to compute the mean
        window spacing Δτ so eigenvalues can be converted into
        continuous-time rates. If omitted, unit spacing is assumed.
    smooth_sigma
        Gaussian smoothing (in grid steps) applied to the |b_k(t)|
        activation magnitudes to suppress high-frequency ringing that
        arises when neighbouring windows differ in eigenvector ordering.
    whiten
        Metric choice for the reduced trajectory. ``False`` (default,
        ship behaviour) keeps ``Z = Y V_r = U_r Σ_r`` so directions with
        larger singular values dominate the fit. ``True`` orthonormalises
        the trajectory (``Z = U_r``) so the fit is balanced across
        modes. Under exact least squares the two operators are similar
        (same eigenvalues, different eigenvectors); with ``ridge > 0``
        the regulariser acts differently in the two metrics, so
        eigenvalues shift slightly as well. Every geometry descriptor
        (Henrici, reactivity, transient gain, eigvec conditioning)
        depends on this choice — document which was used.
    geometry_n_max
        Horizon ``N`` for the transient-gain sweep
        ``max_{n ∈ [1, N]} ||K^n||_2``. Only used with
        ``return_diagnostics=True``; 30 is enough for r ≤ 30 without
        being noticeably expensive.
    return_diagnostics
        If True, additionally return a diagnostics dict (see below).

    Returns
    -------
    patterns : torch.Tensor of shape (K, D, D)
        Real-valued operator archetypes obtained by reshaping
        ``P_r v_k`` back to D x D. Complex conjugate pairs contribute two
        adjacent real matrices (Re, Im) — this preserves phase without
        introducing complex arithmetic to downstream code.
    activations : torch.Tensor of shape (T, K)
        Non-negative mode magnitudes |b_k(t)| where b_k(t) is the
        projection of the reduced trajectory onto the left eigenvector
        of the (local or global) Koopman operator at t. Smoothed by
        ``smooth_sigma``.
    error : float
        Frobenius reconstruction error of the reduced trajectory by the
        top-K modes (average of per-window errors in local mode).
    diagnostics : dict, only when ``return_diagnostics=True``
        - ``eigenvalues`` (T, K) complex λ_k(τ) matched to reference modes.
        - ``growth_rates`` (T, K) Re(log λ / Δτ). Frequency-agnostic and
          identified only up to Δτ discretisation.
        - ``freqs`` (T, K) |Im(log λ / Δτ)| in radians / unit pseudotime.
          Identified only below the Nyquist limit — see ``nyquist_ang_freq``.
        - ``delta_tau`` (T,) mean pseudotime step in the window centered at t
          (== ``dt_seq``; the pair travels with the frequency descriptors so
          the ambiguity caveat lives with the numbers rather than only in the
          docs).
        - ``nyquist_ang_freq`` (T,) π / Δτ, in radians / unit pseudotime.
          Any reported ``freqs[t, k]`` above this value is *aliased*, not
          measured. Compare against ``freqs`` before interpreting oscillations.
        - ``nyquist_cycle_freq`` (T,) 1 / (2 Δτ), same limit in cycles / unit
          pseudotime (convenient when comparing to physical periods).
        - ``geometry`` per-window matrix-geometry descriptors, all
          branch-free (do not depend on log λ):
            * ``henrici``        (T,) departure-from-normality index in [0, 1]
            * ``reactivity``     (T,) σ_max(K) / ρ(K), ≥ 1
            * ``transient_gain`` (T,) max_n ||K^n||_2 over n ≤ ``geometry_n_max``
            * ``eigvec_cond``    (T,) condition number of the eigenvector matrix
          For ``mode='global'`` every entry along T is identical (same K).
        - ``V_r`` (D², r) SVD basis used for the reduction.
        - ``singular_values`` (r,) singular values of the SVD reduction.
        - ``K_seq`` list of (r, r) operators (local mode) or a single
          (r, r) operator (global mode).
        - ``K_ref`` reference operator used to define the pattern basis
          (global K in local mode, or K_seq[0] in global mode).
        - ``dt_seq`` alias of ``delta_tau`` kept for backward compatibility.
        - ``ref_eigenvalues`` (K,) discrete-time eigenvalues of ``K_ref``,
          |λ|-sorted; splitting into (Re, Im) real modes matches ``patterns``.
        - ``mode`` echo of the mode argument.
        - ``reduced_rank`` r used for the SVD reduction.
        - ``whiten`` echo — records the metric under which every
          geometry descriptor above was computed.
    """
    if J_tensor.dim() == 4:
        J_tensor = J_tensor.mean(dim=1)
    if J_tensor.dim() != 3:
        raise ValueError(f"J_tensor must be (T, D, D) or (T, N, D, D); got {tuple(J_tensor.shape)}")

    T, d1, d2 = J_tensor.shape
    if T < 3:
        raise ValueError(f"Koopman decomposition needs at least 3 time points, got T={T}")

    D2 = d1 * d2
    Y = J_tensor.reshape(T, D2).cpu().numpy().astype(np.float64)

    r_default = max(2, min(T - 1, 4 * rank + 8, D2))
    r = int(reduced_rank if reduced_rank is not None else r_default)
    r = max(2, min(r, T - 1, D2))

    Z, V_r, S_r, lift = _reduce_trajectory(Y, r, whiten=whiten)

    if grid is None:
        grid_np = np.arange(T, dtype=np.float64)
    else:
        grid_np = np.asarray(grid, dtype=np.float64)
        if grid_np.shape != (T,):
            raise ValueError(f"grid must have shape ({T},), got {grid_np.shape}")

    mode = mode.lower()
    if mode not in {"local", "global"}:
        raise ValueError(f"mode must be 'local' or 'global', got {mode!r}")

    # ── Fit operator(s) ────────────────────────────────────────────────────
    if mode == "local":
        K_seq, dt_seq = _local_koopman(Z, grid_np, window_half, ridge)
    else:
        K_g = _global_koopman(Z, ridge)
        K_seq = [K_g for _ in range(T)]
        dt_seq = np.full(T, float(np.median(np.diff(grid_np))), dtype=np.float64)

    # ── Extract patterns from a reference operator ─────────────────────────
    # Rationale: local Koopman gives one operator per t; for a stable
    # "archetype" set we take the eigenvectors of the *global* operator
    # (which represents the average dynamics) — or the median-time local
    # operator in 'local' mode. Then per-t amplitudes are computed by
    # projecting z(t) onto the corresponding left eigenvectors. This
    # keeps patterns temporally-stable and activations temporally-varying
    # in the same spirit as semi-NMF's (patterns × activations)
    # factorisation.
    if mode == "local":
        K_ref = _global_koopman(Z, ridge)
        dt_ref = float(np.median(dt_seq[np.isfinite(dt_seq) & (dt_seq > 0)])) \
            if np.any(np.isfinite(dt_seq) & (dt_seq > 0)) else 1.0
    else:
        K_ref = K_seq[0]
        dt_ref = float(dt_seq[0]) if dt_seq[0] > 0 else 1.0

    eigvals_ref, eigvecs_ref = np.linalg.eig(K_ref)
    lam_all, V_real, _ = _split_complex_pairs(eigvals_ref, eigvecs_ref)

    K_out = min(rank, V_real.shape[1])
    if K_out < rank:
        warnings.warn(
            f"Koopman produced only {K_out} real modes (requested rank={rank}). "
            "The reference operator has fewer non-degenerate directions; "
            "increase reduced_rank or check that T is large enough.",
            UserWarning, stacklevel=2,
        )
    lam_top = lam_all[:K_out]
    V_top = V_real[:, :K_out]                                       # (r, K)

    # Lift the reduced eigenvectors back to D² and reshape. ``lift`` is
    # V_r when whiten=False and V_r Σ_r when whiten=True, so the same
    # code path is correct under either metric.
    patterns_flat = lift @ V_top                                    # (D², K)
    patterns_np = patterns_flat.T.reshape(K_out, d1, d2)            # (K, D, D)
    # Unit-norm per pattern so activations carry the scale (mirrors
    # _rescale in semi-NMF and keeps magnitudes comparable across runs).
    p_norm = np.linalg.norm(patterns_np.reshape(K_out, -1), axis=1)
    p_norm = np.where(p_norm < 1e-12, 1.0, p_norm)
    patterns_np = patterns_np / p_norm[:, None, None]
    V_top_scaled = V_top * p_norm[None, :]                          # b_k inherits ||pattern||

    # ── Amplitudes b_k(t) via left-eigenvector projection ─────────────────
    # Compute left eigenvectors of K_ref: solve K_ref^T W = W diag(λ),
    # then normalise so W_k^H V_top_k = 1. This gives biorthogonal
    # projection onto each Koopman mode.
    left_eigvals, W = np.linalg.eig(K_ref.T)
    # Match left eigenvectors to the (already sorted) V_top by nearest
    # eigenvalue on the discrete spectrum.
    W_matched = np.zeros((V_top.shape[0], K_out), dtype=np.complex128)
    used = np.zeros(len(left_eigvals), dtype=bool)
    for k in range(K_out):
        target = lam_top[k]
        # Prefer the exact conjugate structure — match by eigenvalue distance.
        dists = np.where(used, np.inf, np.abs(left_eigvals - target))
        j = int(np.argmin(dists))
        used[j] = True
        w = W[:, j]
        pair = V_top_scaled[:, k].astype(np.complex128)
        norm = np.vdot(w, pair)
        if abs(norm) > 1e-12:
            w = w / norm
        W_matched[:, k] = w

    # b_k(t) = w_k^H z(t) (complex) → take magnitude for non-negativity.
    B = (Z.astype(np.complex128) @ np.conj(W_matched))              # (T, K)
    activations_np = np.abs(B).astype(np.float32)

    if smooth_sigma and smooth_sigma > 0:
        from scipy.ndimage import gaussian_filter1d
        activations_np = gaussian_filter1d(activations_np, sigma=float(smooth_sigma), axis=0)

    # ── Reconstruction error on the reduced trajectory ────────────────────
    Z_hat = np.real(B @ V_top_scaled.T.astype(np.complex128))       # back to r-D
    err = float(np.linalg.norm(Z - Z_hat, ord="fro"))

    patterns = torch.from_numpy(patterns_np.astype(np.float32))
    activations = torch.from_numpy(activations_np)

    if not return_diagnostics:
        return patterns, activations, err

    # ── Local eigenvalues and continuous-time rates per grid point ────────
    T_grid = T
    eig_local = np.full((T_grid, K_out), np.nan, dtype=np.complex128)
    growth = np.full((T_grid, K_out), np.nan, dtype=np.float64)
    freq = np.full((T_grid, K_out), np.nan, dtype=np.float64)
    for t in range(T_grid):
        lam_t, _ = np.linalg.eig(K_seq[t])
        # Match to reference modes by nearest-eigenvalue (Hungarian-lite:
        # greedy nearest neighbour is stable for well-separated spectra).
        used_t = np.zeros(len(lam_t), dtype=bool)
        dt_t = dt_seq[t] if dt_seq[t] > 0 else dt_ref
        for k in range(K_out):
            dists = np.where(used_t, np.inf, np.abs(lam_t - lam_top[k]))
            j = int(np.argmin(dists))
            used_t[j] = True
            eig_local[t, k] = lam_t[j]
            # log of a complex eigenvalue; guard against |λ| ≈ 0.
            if abs(lam_t[j]) > 1e-12 and dt_t > 0:
                mu = np.log(lam_t[j]) / dt_t
                growth[t, k] = float(mu.real)
                freq[t, k] = float(abs(mu.imag))

    # ── Nyquist limits (per window, in the units of ``freqs``) ────────────
    # freqs is stored in rad / unit-τ, so the matching Nyquist bound is
    # π / Δτ. The cycles-per-unit-τ variant is exposed too, since the
    # user often prefers to reason in physical periods.
    safe_dt = np.where(dt_seq > 0, dt_seq, np.nan)
    nyquist_ang = np.pi / safe_dt
    nyquist_cyc = 1.0 / (2.0 * safe_dt)

    # ── Per-window branch-free geometry descriptors ───────────────────────
    henrici_arr = np.zeros(T_grid, dtype=np.float64)
    reactivity_arr = np.zeros(T_grid, dtype=np.float64)
    trans_gain_arr = np.zeros(T_grid, dtype=np.float64)
    eigcond_arr = np.zeros(T_grid, dtype=np.float64)
    if mode == "global":
        g = _matrix_geometry(K_seq[0], n_max=geometry_n_max)
        henrici_arr[:] = g["henrici"]
        reactivity_arr[:] = g["reactivity"]
        trans_gain_arr[:] = g["transient_gain"]
        eigcond_arr[:] = g["eigvec_cond"]
    else:
        for t in range(T_grid):
            g = _matrix_geometry(K_seq[t], n_max=geometry_n_max)
            henrici_arr[t] = g["henrici"]
            reactivity_arr[t] = g["reactivity"]
            trans_gain_arr[t] = g["transient_gain"]
            eigcond_arr[t] = g["eigvec_cond"]

    diagnostics = {
        "eigenvalues": eig_local,
        "growth_rates": growth,
        "freqs": freq,
        "delta_tau": dt_seq,
        "nyquist_ang_freq": nyquist_ang,
        "nyquist_cycle_freq": nyquist_cyc,
        "geometry": {
            "henrici":        henrici_arr,
            "reactivity":     reactivity_arr,
            "transient_gain": trans_gain_arr,
            "eigvec_cond":    eigcond_arr,
        },
        "V_r": V_r,
        "singular_values": S_r,
        "K_seq": K_seq if mode == "local" else K_seq[0],
        "K_ref": K_ref,
        "dt_seq": dt_seq,                # alias, kept for backward compat
        "mode": mode,
        "reduced_rank": r,
        "window_half": window_half if mode == "local" else None,
        "ref_eigenvalues": lam_top,
        "whiten": bool(whiten),
    }
    return patterns, activations, err, diagnostics
