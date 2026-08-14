"""
High-level drift field training and Jacobian tensor decomposition.
All results are stored in adata.uns[key_added].
"""
from __future__ import annotations

import warnings
from typing import Optional

import numpy as np
import torch
from scipy.ndimage import gaussian_filter1d
from tqdm.auto import trange


# ---------------------------------------------------------------------------
# Pseudotime velocity helper
# ---------------------------------------------------------------------------

def _pseudotime_velocity(X_pca: np.ndarray, pseudotime: np.ndarray, k: int = 15) -> np.ndarray:
    """
    Compute a directional velocity from DPT pseudotime.
    For each cell: weighted average displacement toward higher-pseudotime neighbours.
    Gives the drift field a directional signal when RNA velocity is unavailable.
    """
    from sklearn.neighbors import NearestNeighbors
    N = X_pca.shape[0]
    V = np.zeros_like(X_pca, dtype=np.float32)
    nbrs = NearestNeighbors(n_neighbors=k + 1, metric="euclidean").fit(X_pca)
    _, idx = nbrs.kneighbors(X_pca)
    for i in range(N):
        nbr = idx[i, 1:]
        dt  = pseudotime[nbr] - pseudotime[i]
        fwd = dt > 0
        if not fwd.any():
            continue
        w    = dt[fwd]
        disp = X_pca[nbr[fwd]] - X_pca[i]
        V[i] = (w[:, None] * disp).sum(0) / w.sum()
    mag = np.linalg.norm(V, axis=1).mean() + 1e-8
    return V / mag   # unit-mean magnitude


def _biased_velocity(
    X_pca: np.ndarray,
    pseudotime: np.ndarray,
    terminal_centroid: np.ndarray,
    progenitor_centroid: Optional[np.ndarray] = None,
    bias_strength: float = 1.0,
    k: int = 15,
) -> np.ndarray:
    """
    Compute a velocity prior biased toward a known terminal cell state.

    Combines two components:

    1. **Pseudotime gradient** — points each cell toward higher-pseudotime
       neighbours (same as ``_pseudotime_velocity``).

    2. **Target pull** — points each cell toward the terminal state centroid,
       scaled by pseudotime so committed cells are pulled harder than progenitors:

           pull[i] = normalize(X_terminal - X[i]) × pseudotime[i]

    Result:  V = V_pseudotime  +  bias_strength × target_pull

    The combined vector is then normalised to unit-mean magnitude before being
    passed to the drift field as the velocity prior.

    Parameters
    ----------
    X_pca              : (N, D) cell coordinates in latent space.
    pseudotime         : (N,)  pseudotime values in [0, 1].
    terminal_centroid  : (D,)  centroid of the terminal cell population.
    progenitor_centroid: (D,)  optional centroid of the progenitor population.
                               If provided, the pull is measured relative to
                               the progenitor-to-terminal axis.
    bias_strength      : Scaling factor for the target pull relative to the
                         pseudotime gradient.
                         ``0``  = pure pseudotime gradient (no bias).
                         ``1``  = equal weight (default, recommended starting point).
                         ``2-5``= strong bias toward terminal (use when the
                                  pseudotime gradient is noisy at the branch point).
    k                  : Neighbours for pseudotime gradient computation.

    Examples
    --------
    >>> X_ery    = ad.obsm['X_pca'][ery_mask]
    >>> pt_ery   = ad.obs['pseudotime'].values[ery_mask]
    >>> X_term   = ad.obsm['X_pca'][ad.obs['paul15_clusters']=='3Ery'].mean(0)
    >>> V = _biased_velocity(X_ery, pt_ery, terminal_centroid=X_term,
    ...                      bias_strength=1.5)
    """
    # ── Component 1: pseudotime gradient ──────────────────────────────────
    V_pt = _pseudotime_velocity(X_pca, pseudotime, k=k)  # already unit-mean

    # ── Component 2: target pull ───────────────────────────────────────────
    # Direction from each cell toward terminal centroid
    direction = terminal_centroid[np.newaxis, :] - X_pca   # (N, D)

    # If progenitor centroid is provided, project direction onto the
    # progenitor→terminal axis (removes lateral drift perpendicular to the path)
    if progenitor_centroid is not None:
        axis = terminal_centroid - progenitor_centroid          # (D,)
        axis = axis / (np.linalg.norm(axis) + 1e-8)
        proj = (direction * axis[np.newaxis, :]).sum(1, keepdims=True)  # scalar projection
        direction = proj * axis[np.newaxis, :]                  # projected direction

    # Normalise per cell so direction magnitude doesn't depend on distance
    norms     = np.linalg.norm(direction, axis=1, keepdims=True) + 1e-8
    direction = direction / norms                              # unit vectors

    # Scale by pseudotime: committed cells (high pt) get stronger pull
    # Progenitors (pt≈0) are barely pulled — they follow the gradient instead
    pull = direction * pseudotime[:, np.newaxis]               # (N, D)

    # ── Combine ────────────────────────────────────────────────────────────
    V_combined = V_pt + bias_strength * pull.astype(np.float32)
    mag = np.linalg.norm(V_combined, axis=1).mean() + 1e-8
    return (V_combined / mag).astype(np.float32)


# ---------------------------------------------------------------------------
# Pseudotime window construction
# ---------------------------------------------------------------------------

def _serialize_sweep(sweep) -> Optional[list]:
    """Strip ndarray fields from the bandwidth-sweep table for .uns storage."""
    if sweep is None:
        return None
    out = []
    for row in sweep:
        out.append({
            "h": float(row.get("h", float("nan"))),
            "R": float(row.get("R", float("nan"))),
            "C": float(row.get("C", float("nan"))),
            "L": float(row.get("L", float("nan"))),
            "S": float(row.get("S", float("nan"))),
            "n_eff_interior_min": float(row.get("n_eff_interior_min", float("nan"))),
        })
    return out


def _compute_jacobians_per_cell(model, X: torch.Tensor, T: torch.Tensor,
                                 batch_size: int = 128, verbose: bool = True) -> np.ndarray:
    """Compute model Jacobians for every cell in batches. Returns (N, D, D) numpy."""
    N, D = X.shape
    mem_gb = N * D * D * 4 / 1e9
    if mem_gb > 4.0:
        warnings.warn(
            f"Per-cell Jacobian tensor will use ~{mem_gb:.1f} GB. "
            "Consider windowing='fixed' for very large datasets.",
            ResourceWarning, stacklevel=2,
        )
    Js: list[np.ndarray] = []
    iterator = trange(0, N, batch_size, desc="Per-cell Jacobians",
                      leave=False) if verbose else range(0, N, batch_size)
    with torch.enable_grad():
        for i in iterator:
            xb = X[i:i + batch_size]; tb = T[i:i + batch_size]
            Jb = model.jacobian(xb, tb)
            Js.append(Jb.detach().cpu().numpy().astype(np.float32))
    return np.concatenate(Js, axis=0)


def _build_windows(
    pseudotime: np.ndarray,
    n_windows: int = 100,
    overlap: float = 0.80,
    min_cells: int = 5,
) -> list[tuple[float, float, np.ndarray]]:
    """Build overlapping pseudotime windows with guaranteed minimum cell count."""
    width  = 1.0 / (1.0 + (n_windows - 1) * (1.0 - overlap))
    step   = width * (1.0 - overlap)
    windows = []
    for i in range(n_windows):
        lo  = i * step
        hi  = lo + width
        idx = np.where((pseudotime >= lo) & (pseudotime <= hi))[0]
        # Expand window if too sparse
        if len(idx) < min_cells:
            dists = np.abs(pseudotime - (lo + hi) / 2)
            idx   = np.argsort(dists)[:min_cells]
        windows.append((lo, hi, idx))
    return windows


# ---------------------------------------------------------------------------
# Gene-space projection
# ---------------------------------------------------------------------------

def _gene_scores(patterns: torch.Tensor, pca_loadings: np.ndarray, sign_flip: list[float]):
    """Project archetype patterns to gene space and return ranked gene scores."""
    K = patterns.shape[0]
    results = {}
    for k in range(K):
        A_k     = patterns[k].numpy() * sign_flip[k]
        eigvals, eigvecs = np.linalg.eig(A_k)
        max_idx = int(np.argmax(np.real(eigvals)))
        v_max   = np.real(eigvecs[:, max_idx])
        scores  = pca_loadings @ v_max        # (n_genes,)
        results[k] = scores
    return results


# ---------------------------------------------------------------------------
# Main API function
# ---------------------------------------------------------------------------

def fit_drift(
    adata,
    *,
    time_key: str = "pseudotime",
    rep: Optional[str] = None,          # None = auto-detect from prepare_trajectory
    branch_key: Optional[str] = None,   # obsm key for branch probabilities (Palantir)
    # Model
    hidden: int = 256,
    depth: int = 4,
    beta: float = 0.1,
    use_spectral_norm: bool = True,
    # Velocity prior
    vel_scale: float = 2.0,
    vel_time_mode: str = "flat",
    # Training
    n_epochs: int = 5000,
    batch_size: int = 512,
    lr: float = 2e-4,
    weight_decay: float = 1e-4,
    sigma: float = 0.1,
    use_local_sigma: bool = False,
    seed: int = 42,
    device: Optional[str] = None,
    # Temporal aggregation — default: adaptive Gaussian kernel
    windowing: str = "kernel",                       # "kernel" | "fixed"
    bandwidth = "auto",                              # float | "auto"
    bandwidth_grid: tuple = (0.01, 0.02, 0.03, 0.05, 0.08, 0.10),
    n_eff_min: float = 30.0,
    adaptive: bool = False,                          # kNN-adaptive h(tau)
    knn_k: int = 80,
    grid_size: int = 200,                            # kernel-mode eval grid
    n_boot: int = 20,                                # bootstraps for bandwidth scoring
    # Legacy fixed-window scheme
    n_windows: int = 100,
    overlap: float = 0.80,
    smooth_sigma: float = 1.5,
    # Archetype decomposition
    n_archetypes: int = 5,
    n_restarts: int = 5,
    tv_lambda: float = 0.0,
    archetype_method: str = "snmf",              # "snmf" | "koopman"
    koopman_kwargs: Optional[dict] = None,
    # Output
    key_added: str = "scjdo",
    verbose: bool = True,
):
    """
    Train a drift field and compute archetype decomposition.

    Automatically derives a pseudotime-gradient velocity prior when no RNA
    velocity is available, so the drift field follows the developmental path.

    Results are stored in ``adata.uns[key_added]``:

    ========================  ================================================
    Key                       Content
    ========================  ================================================
    ``J_tensor``              Jacobian tensor (T_eval, D, D), numpy
    ``t_centers``             Pseudotime centers (T_eval,), numpy
    ``patterns``              Archetype patterns (K, D, D), numpy
    ``activations``           Temporal activations (T_eval, K), numpy
    ``max_real_eig``          Max real eigenvalue per time point (T_eval,)
    ``gene_scores``           Gene loadings per archetype {k: (n_genes,)}
    ``corr_mat``              Pairwise activation correlations (K, K)
    ``windowing``             "kernel" or "fixed"
    ``bandwidth``             Selected bandwidth (kernel mode only)
    ``n_eff``                 Effective sample size per grid point (kernel only)
    ``kernel_score``          {R, C, L, S} for the selected bandwidth (kernel)
    ``kernel_sweep``          Per-bandwidth scores (kernel + auto only)
    ``lam_bootstrap``         (n_boot, T_eval) bootstrap λ_max curves (kernel)
    ``params``                Training configuration dict
    ========================  ================================================

    Parameters
    ----------
    adata         : AnnData after ``sjd.pp.prepare_trajectory``.
    time_key      : Column in ``adata.obs`` with normalized pseudotime [0,1].
    rep           : Key in ``adata.obsm`` for the latent representation.
    n_epochs      : Training iterations.
    n_archetypes  : Number of operator archetypes (K).
    vel_scale     : Strength of the pseudotime-gradient velocity prior.
    vel_time_mode : Gate shape — 'flat' (constant), 'root', 'rise', 'mid'.
    windowing     : 'kernel' (default) uses adaptive Gaussian kernel windowing
                    with bootstrap-selected bandwidth. 'fixed' uses the legacy
                    100-window / 80%-overlap scheme.
    bandwidth     : 'auto' (default) selects h* by maximising S=R*C*L subject
                    to n_eff>=n_eff_min. Pass a float to set h directly.
    bandwidth_grid: Candidate h values for the auto sweep.
    n_eff_min     : Floor on effective sample size (cells) per grid point.
    adaptive      : If True, use kNN-adaptive h(τ) instead of a global h.
    knn_k         : Neighbour count for the adaptive bandwidth.
    grid_size     : Number of evaluation points in pseudotime (kernel mode).
    n_boot        : Bootstrap replicates for reproducibility scoring.
    n_windows     : Number of fixed windows (used only when windowing='fixed').
    overlap       : Fixed-window overlap fraction (legacy).
    smooth_sigma  : Post-binning Gaussian smoothing sigma (legacy).
    archetype_method : Decomposition backend for the temporal Jacobian tensor.
                    ``'snmf'`` (default) uses non-negative activations of signed
                    pattern operators (interpretable "programs"). ``'koopman'``
                    uses windowed EDMD on the vectorised trajectory: modes are
                    Koopman eigenvectors and per-time amplitudes come from the
                    left-eigenvector projection, giving growth rates and
                    oscillation frequencies as spectral summaries of the
                    operator sequence. Best for dense, time-resolved,
                    mostly-autonomous trajectories; less suited to sparse,
                    branch-heavy atlases. See
                    :func:`scjdo.archetypes.koopman.koopman_modes`.
    koopman_kwargs : Extra keyword arguments forwarded to
                    :func:`koopman_modes` when ``archetype_method='koopman'``.
                    Common overrides: ``reduced_rank``, ``mode='local'|'global'``,
                    ``window_half``, ``ridge``, ``smooth_sigma``.
    key_added     : Where to store results in ``adata.uns``. Koopman-specific
                    spectral diagnostics (eigenvalues, growth rates, oscillation
                    frequencies, per-window operators) are stored under
                    ``adata.uns[key_added]['koopman']`` when
                    ``archetype_method='koopman'``.

    Examples
    --------
    >>> sjd.tl.fit_drift(adata, time_key="pseudotime", n_archetypes=5, n_epochs=5000)
    >>> adata.uns["scjdo"].keys()
    dict_keys(['J_tensor', 't_centers', 'patterns', 'activations', ...])
    """
    from scjdo.models.drift import DriftField, DriftConfig
    from scjdo.losses import denoising_score_matching, control_energy, local_sigma
    from scjdo.archetypes.decompose import jacobian_modes

    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"

    torch.manual_seed(seed)
    np.random.seed(seed)

    # ── Auto-detect representation from prepare_trajectory metadata ───────
    if rep is None:
        prep = adata.uns.get("scjdo_prep", {})
        rep  = prep.get("rep", "X_pca")
        if verbose:
            print(f"[fit_drift] Using representation: {rep} "
                  f"(latent={prep.get('latent','pca')})")

    # ── Extract tensors ────────────────────────────────────────────────────
    if rep not in adata.obsm:
        raise KeyError(f"Representation '{rep}' not in adata.obsm. Run sjd.pp.prepare_trajectory first.")
    if time_key not in adata.obs.columns:
        raise KeyError(f"time_key='{time_key}' not in adata.obs. Run sjd.pp.prepare_trajectory first.")

    X_np = adata.obsm[rep].astype(np.float32)
    T_np = adata.obs[time_key].values.astype(np.float32)
    N, D = X_np.shape

    X = torch.tensor(X_np, device=device)
    T = torch.tensor(T_np, device=device)

    # ── Pseudotime-gradient velocity (auto) ───────────────────────────────
    if verbose:
        print("Computing pseudotime-gradient velocity prior...")
    V_np     = _pseudotime_velocity(X_np, T_np, k=15)
    V_tensor = torch.tensor(V_np, device=device)

    # ── Build model ────────────────────────────────────────────────────────
    cfg = DriftConfig(
        dim=D, hidden=hidden, depth=depth, beta=beta,
        use_spectral_norm=use_spectral_norm,
        use_velocity_prior=True,
        vel_scale=vel_scale,
        vel_k=15,
        vel_time_mode=vel_time_mode,
    )
    model = DriftField(cfg, X_ref=X, V_ref=V_tensor).to(device)
    opt   = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=n_epochs)

    if verbose:
        print(f"DriftField: {sum(p.numel() for p in model.parameters()):,} parameters | "
              f"device={device} | epochs={n_epochs}")

    # Precompute local sigma if requested
    sigma_per_cell = None
    if use_local_sigma:
        sigma_per_cell = local_sigma(X).cpu()

    # ── Training ───────────────────────────────────────────────────────────
    losses = []
    model.train()
    for step in trange(n_epochs, desc="Training drift field", disable=not verbose):
        idx    = torch.randint(0, N, (batch_size,), device=device)
        xb, tb = X[idx], T[idx]

        if sigma_per_cell is not None:
            sig = sigma_per_cell[idx.cpu()].to(device)
        else:
            sig = sigma

        loss = denoising_score_matching(model, xb, tb, sigma=sig)
        loss = loss + cfg.alpha_control * control_energy(model(xb, tb))

        opt.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step()
        sched.step()
        losses.append(loss.item())

    model.eval()
    if verbose:
        print(f"Training complete. Final loss: {losses[-1]:.4f}")

    # ── Branch-probability weights (Palantir or similar) ──────────────────
    branch_weights = None
    if branch_key is not None:
        if branch_key in adata.obsm:
            bp = adata.obsm[branch_key].astype(np.float32)
            branch_weights = bp[:, np.argmax(bp.mean(0))]
            if verbose:
                branch_names = adata.uns.get("palantir_branch_names",
                                             [f"branch_{i}" for i in range(bp.shape[1])])
                dominant = branch_names[np.argmax(bp.mean(0))]
                print(f"[fit_drift] Branch weighting: dominant branch = {dominant}")
        elif branch_key in adata.obs.columns:
            branch_weights = adata.obs[branch_key].values.astype(np.float32)
            if verbose:
                print(f"[fit_drift] Branch weighting from adata.obs['{branch_key}']")
        else:
            warnings.warn(
                f"branch_key='{branch_key}' not found in adata.obsm or adata.obs. "
                "Running without branch weighting.", UserWarning, stacklevel=2
            )

    # ── Temporal aggregation ───────────────────────────────────────────────
    windowing = str(windowing).lower()
    kernel_diag: dict = {}

    if windowing == "kernel":
        if verbose:
            print("Computing per-cell Jacobians...")
        J_per_cell = _compute_jacobians_per_cell(
            model, X, T, batch_size=128, verbose=verbose,
        )
        if verbose:
            print(f"Kernel-weighted temporal operator: grid_size={grid_size}, "
                  f"bandwidth={bandwidth!r}, adaptive={adaptive}")
        from scjdo.archetypes.windowing import build_temporal_operator
        kw = build_temporal_operator(
            J_per_cell, T_np,
            grid_size=grid_size,
            bandwidth=bandwidth,
            bandwidth_grid=tuple(bandwidth_grid),
            n_eff_min=n_eff_min,
            n_boot=n_boot,
            adaptive=adaptive,
            knn_k=knn_k,
            seed=seed,
            cell_weights=branch_weights,
            verbose=verbose,
        )
        J_np      = kw["J_tensor"]
        t_centers = kw["grid"]
        kernel_diag = {
            "bandwidth":  kw["bandwidth"],
            "n_eff":      kw["n_eff"],
            "lam_max":    kw["lam"],
            "score":      kw["score"],
            "sweep":      kw["sweep"],
            "boots":      kw["boots"],
            "mode":       kw["mode"],
        }
        # Free per-cell tensor early
        del J_per_cell

    elif windowing == "fixed":
        if verbose:
            print(f"Fixed-window temporal aggregation: n_windows={n_windows}, "
                  f"overlap={overlap}, smooth_sigma={smooth_sigma}")
        windows   = _build_windows(T_np, n_windows=n_windows, overlap=overlap)
        J_tensor  = torch.zeros(n_windows, D, D)
        t_centers = np.zeros(n_windows, dtype=np.float32)
        iter_j = trange(n_windows, desc="Jacobians", leave=False) if verbose else range(n_windows)
        for i in iter_j:
            lo, hi, idx_w = windows[i]
            t_centers[i] = (lo + hi) / 2.0
            xb = X[idx_w].to(device)
            tb = T[idx_w].to(device)
            with torch.enable_grad():
                Jb = model.jacobian(xb, tb)
            if branch_weights is not None:
                w = torch.tensor(branch_weights[idx_w], device=Jb.device)
                w = w / (w.sum() + 1e-8)
                J_tensor[i] = (Jb.detach().cpu() * w.cpu().view(-1, 1, 1)).sum(0)
            else:
                J_tensor[i] = Jb.detach().cpu().mean(dim=0)
        J_np = gaussian_filter1d(J_tensor.numpy(), sigma=smooth_sigma, axis=0)

    else:
        raise ValueError(f"windowing must be 'kernel' or 'fixed', got {windowing!r}")

    T_eval = J_np.shape[0]

    archetype_method = str(archetype_method).lower()
    if archetype_method not in {"snmf", "koopman"}:
        raise ValueError(f"archetype_method must be 'snmf' or 'koopman', got {archetype_method!r}")

    if verbose:
        print(f"Running archetype decomposition ({archetype_method})...")

    # ── Archetype decomposition ────────────────────────────────────────────
    J_smooth = torch.from_numpy(J_np)
    koopman_diag: Optional[dict] = None
    if archetype_method == "snmf":
        patterns, activations, recon_err = jacobian_modes(
            J_smooth, rank=n_archetypes, n_restarts=n_restarts, seed=seed,
            tv_lambda=tv_lambda,
        )
    else:
        from scjdo.archetypes.koopman import koopman_modes
        kkw = dict(koopman_kwargs or {})
        kkw.setdefault("mode", "local")
        # Pass the pseudotime grid so eigenvalues convert to continuous-time rates.
        patterns, activations, recon_err, koopman_diag = koopman_modes(
            J_smooth, rank=n_archetypes,
            grid=t_centers,
            return_diagnostics=True,
            **kkw,
        )
    activations_np = activations.numpy()   # (T, K) — non-negative

    # The Koopman backend may return fewer real modes than requested
    # (complex conjugate splitting, degenerate spectra). Use the actual
    # returned count downstream so shape assumptions stay consistent.
    n_archetypes_eff = int(patterns.shape[0])

    # Variance explained (R²)
    total_sq = float((J_smooth.reshape(T_eval, -1) ** 2).sum())
    recon_sq = float(np.sum((J_np.reshape(T_eval, -1) -
                             (activations_np @ patterns.numpy().reshape(n_archetypes_eff, -1))) ** 2))
    r2 = max(0.0, 1.0 - recon_sq / (total_sq + 1e-8))

    # Max real eigenvalue per time point (local sensitivity)
    if kernel_diag.get("lam_max") is not None:
        max_eig = np.asarray(kernel_diag["lam_max"], dtype=np.float32)
    else:
        max_eig = np.array([
            float(np.real(np.linalg.eigvals(J_np[i])).max())
            for i in range(T_eval)
        ])

    # Normalize activations for downstream plotting
    act_norm = activations_np.copy()
    for k in range(n_archetypes_eff):
        a = act_norm[:, k]
        act_norm[:, k] = a / (a.max() + 1e-8)

    # Sign convention: activations are already ≥ 0 from semi-NMF / Koopman |b_k|.
    sign_flip = [1.0] * n_archetypes_eff

    # Pairwise temporal correlations
    corr_mat = np.corrcoef(act_norm.T)   # (K, K)

    # ── Gene-space projection (per archetype) ─────────────────────────────
    gene_scores_dict = {}
    top_genes_dict   = {}
    pca_loadings     = None
    gene_names       = []
    if "PCs" in adata.varm:
        pca_loadings = adata.varm["PCs"].astype(np.float32)   # (n_genes, n_pcs)
        gene_names   = list(adata.var_names)
        gene_scores_dict = _gene_scores(patterns, pca_loadings, sign_flip)
        for k, gs in gene_scores_dict.items():
            ranked = np.argsort(np.abs(gs))[::-1]
            top_genes_dict[k] = [gene_names[j] for j in ranked[:50]]
            gene_scores_dict[k] = gs

    # ── Per-window instability gene scores ─────────────────────────────────
    # For each sensitive window (Re(λ_max) > 0), project the eigenvector
    # associated with max Re(λ) to gene space. Sign consistency is enforced
    # by aligning each window's eigenvector to the previous one.
    instability_scores  = np.zeros((T_eval, len(gene_names)), dtype=np.float32)
    instability_evecs   = np.zeros((T_eval, D), dtype=np.float32)
    prev_v = None
    for i in range(T_eval):
        if max_eig[i] <= 0 or pca_loadings is None:
            prev_v = None
            continue
        eigvals, eigvecs = np.linalg.eig(J_np[i])
        max_idx  = int(np.argmax(np.real(eigvals)))
        v        = np.real(eigvecs[:, max_idx])
        v        = v / (np.linalg.norm(v) + 1e-8)
        # Flip sign to be consistent with previous sensitive window
        if prev_v is not None and np.dot(v, prev_v) < 0:
            v = -v
        prev_v = v
        instability_evecs[i]  = v
        instability_scores[i] = pca_loadings @ v          # (n_genes,)

    # Global instability ranking: mean |score| across sensitive windows
    sensitive_mask = max_eig > 0
    if sensitive_mask.sum() > 0 and len(gene_names) > 0:
        mean_instab = np.abs(instability_scores[sensitive_mask]).mean(0)
        top_instab_idx   = np.argsort(mean_instab)[::-1]
        top_instab_genes = [gene_names[j] for j in top_instab_idx[:50]]
        top_instab_scores = mean_instab[top_instab_idx[:50]]
    else:
        top_instab_genes  = []
        top_instab_scores = np.array([])

    # Per-archetype instability genes: which genes drive instability
    # specifically when archetype k is dominant?
    arch_instab_genes = {}
    arch_instab_scores = {}
    for k in range(n_archetypes_eff):
        # Windows where this archetype is in top quartile of activation
        thresh    = np.quantile(act_norm[:, k], 0.75)
        arch_mask = (act_norm[:, k] >= thresh) & sensitive_mask
        if arch_mask.sum() > 0 and len(gene_names) > 0:
            arch_mean = np.abs(instability_scores[arch_mask]).mean(0)
            ranked    = np.argsort(arch_mean)[::-1]
            arch_instab_genes[str(k)]  = [gene_names[j] for j in ranked[:50]]
            arch_instab_scores[str(k)] = arch_mean[ranked[:50]]
        else:
            arch_instab_genes[str(k)]  = []
            arch_instab_scores[str(k)] = np.array([])

    # ── Store results ──────────────────────────────────────────────────────
    adata.uns[key_added] = {
        "J_tensor":    J_np,
        "t_centers":   t_centers,
        "patterns":    patterns.numpy(),
        "activations": activations_np,
        "act_norm":    act_norm,
        "max_real_eig": max_eig,
        "corr_mat":    corr_mat,
        "gene_scores":  {str(k): v for k, v in gene_scores_dict.items()},
        "top_genes":    {str(k): v for k, v in top_genes_dict.items()},
        # instability-specific outputs
        "instability_scores":      instability_scores,    # (n_windows, n_genes)
        "top_instability_genes":   top_instab_genes,      # list[str], global ranking
        "top_instability_scores":  top_instab_scores,     # mean |score| per gene
        "arch_instability_genes":  arch_instab_genes,     # per-archetype ranking
        "arch_instability_scores": arch_instab_scores,
        "gene_names":  gene_names,
        "losses":      losses,
        "r2":          r2,
        "recon_err":   recon_err,
        "params": {
            "n_archetypes": n_archetypes,
            "n_archetypes_eff": n_archetypes_eff,
            "archetype_method": archetype_method,
            "windowing": windowing,
            "n_windows": n_windows,          # legacy fixed-mode kwarg
            "T_eval": T_eval,                # actual length of J_tensor / t_centers
            "grid_size": grid_size,          # kernel-mode eval grid
            "bandwidth": kernel_diag.get("bandwidth"),
            "n_eff_min": n_eff_min,
            "adaptive": adaptive,
            "knn_k": knn_k,
            "n_boot": n_boot,
            "overlap": overlap,
            "smooth_sigma": smooth_sigma,
            "n_epochs": n_epochs, "hidden": hidden, "depth": depth,
            "beta": beta, "vel_scale": vel_scale, "vel_time_mode": vel_time_mode,
            "lr": lr, "sigma": sigma, "seed": seed, "rep": rep, "time_key": time_key,
        },
        # Kernel-windowing diagnostics (None / empty when windowing=='fixed')
        "windowing":      windowing,
        "bandwidth":      kernel_diag.get("bandwidth"),
        "n_eff":          kernel_diag.get("n_eff"),
        "kernel_score":   kernel_diag.get("score"),
        "kernel_sweep":   _serialize_sweep(kernel_diag.get("sweep")),
        "lam_bootstrap":  kernel_diag.get("boots"),
    }

    # Attach Koopman spectral diagnostics under a sub-key when applicable.
    # These are complementary to the (patterns, activations) that already
    # slot into the semi-NMF-shaped fields above.
    if koopman_diag is not None:
        adata.uns[key_added]["koopman"] = {
            "eigenvalues":        koopman_diag.get("eigenvalues"),
            "growth_rates":       koopman_diag.get("growth_rates"),
            "freqs":              koopman_diag.get("freqs"),
            # Branch-ambiguity companions — freqs above these limits are aliased.
            "delta_tau":          koopman_diag.get("delta_tau"),
            "nyquist_ang_freq":   koopman_diag.get("nyquist_ang_freq"),
            "nyquist_cycle_freq": koopman_diag.get("nyquist_cycle_freq"),
            # Branch-free geometry (does not depend on log λ) — often the
            # better-posed half of the output in the snapshot setting.
            "geometry":           koopman_diag.get("geometry"),
            "singular_values":    koopman_diag.get("singular_values"),
            "dt_seq":             koopman_diag.get("dt_seq"),
            "ref_eigenvalues":    koopman_diag.get("ref_eigenvalues"),
            "reduced_rank":       koopman_diag.get("reduced_rank"),
            "window_half":        koopman_diag.get("window_half"),
            "mode":               koopman_diag.get("mode"),
            "whiten":             koopman_diag.get("whiten"),
        }

    # Store model drift and pseudotime velocity for visualization
    with torch.no_grad():
        drift_all = model(X, T).cpu().numpy()
    adata.obsm["X_drift"]          = drift_all   # full model drift (50D)
    adata.obsm["X_velocity_pseudo"] = V_np        # pseudotime-gradient velocity (50D)

    if verbose:
        print(f"Done. R²={r2:.3f} | Archetypes stored in adata.uns['{key_added}']")

    return model


# ---------------------------------------------------------------------------
# Gene extraction utility
# ---------------------------------------------------------------------------

def get_instability_genes(
    adata,
    key: str = "scjdo",
    n_genes: int = 20,
    activation_threshold: float = 0.5,
    sensitivity_threshold: float = 0.05,
    min_sensitive_fraction: float = 0.0,
    top_archetypes: Optional[int] = None,
) -> "pd.DataFrame":
    """
    Extract the top genes driving instability from each archetype,
    ranked by how much the archetype is associated with locally sensitive
    (positive-eigenvalue) pseudotime windows.

    An archetype is characterised as "instable" when it is strongly active
    (activation > ``activation_threshold``) during windows where the
    maximum real eigenvalue exceeds ``sensitivity_threshold``.

    Parameters
    ----------
    n_genes               : Top N genes to return per archetype.
    activation_threshold  : Minimum normalised activation to consider a
                            window "active" for that archetype.
    sensitivity_threshold : Minimum Re(λ_max) to consider a window sensitive.
    min_sensitive_fraction: Minimum fraction of active windows that must be
                            sensitive to include the archetype (0 = all).
    top_archetypes        : If set, return only the K most instable archetypes.

    Returns
    -------
    pandas.DataFrame with columns:
        archetype, mean_sensitivity, sensitive_fraction,
        peak_pseudotime, gene, instability_score, rank
    """
    import pandas as pd

    res         = adata.uns[key]
    eig         = res["max_real_eig"]          # (n_windows,)
    act_norm    = res["act_norm"]              # (n_windows, K)
    t_np        = res["t_centers"]
    gene_names  = res["gene_names"]
    instab_sc   = res["instability_scores"]    # (n_windows, n_genes)
    arch_genes  = res["arch_instability_genes"]
    arch_scores = res["arch_instability_scores"]

    if not gene_names:
        raise ValueError("No gene names. Re-run sjd.tl.fit_drift on adata with PCA loadings.")

    K = act_norm.shape[1]
    sens_mask = eig > sensitivity_threshold

    # ── Score each archetype by mean sensitivity during active windows ──────
    arch_stats = []
    for k in range(K):
        active_mask = act_norm[:, k] > activation_threshold
        n_active    = int(active_mask.sum())
        if n_active == 0:
            continue
        both_mask       = active_mask & sens_mask
        mean_sens       = float(eig[active_mask].mean())
        sens_frac       = float(both_mask.sum() / n_active)
        peak_pt         = float(t_np[np.argmax(act_norm[:, k])])

        if sens_frac < min_sensitive_fraction:
            continue

        arch_stats.append({
            "k":                  k,
            "archetype":          f"A{k+1}",
            "mean_sensitivity":   round(mean_sens, 4),
            "sensitive_fraction": round(sens_frac, 3),
            "peak_pseudotime":    round(peak_pt, 3),
        })

    if not arch_stats:
        raise ValueError("No archetypes passed the sensitivity filter. "
                         "Lower min_sensitive_fraction or sensitivity_threshold.")

    # Sort by mean sensitivity (most instable first)
    arch_stats.sort(key=lambda x: x["mean_sensitivity"], reverse=True)

    if top_archetypes is not None:
        arch_stats = arch_stats[:top_archetypes]

    # ── Build gene table ───────────────────────────────────────────────────
    rows = []
    for stat in arch_stats:
        k        = stat["k"]
        genes_k  = arch_genes.get(str(k), [])[:n_genes]
        scores_k = arch_scores.get(str(k), np.array([]))[:n_genes]

        for rank, (gene, score) in enumerate(zip(genes_k, scores_k), start=1):
            rows.append({
                "archetype":          stat["archetype"],
                "mean_sensitivity":   stat["mean_sensitivity"],
                "sensitive_fraction": stat["sensitive_fraction"],
                "peak_pseudotime":    stat["peak_pseudotime"],
                "gene":               gene,
                "instability_score":  round(float(score), 4),
                "rank":               rank,
            })

    df = pd.DataFrame(rows, columns=[
        "archetype", "mean_sensitivity", "sensitive_fraction",
        "peak_pseudotime", "gene", "instability_score", "rank",
    ])
    return df


# ---------------------------------------------------------------------------
# Branch-separated analysis
# ---------------------------------------------------------------------------

def fit_drift_branches(
    adata,
    *,
    branch_key: str = "palantir_branch_probs",
    branch_names: Optional[list[str]] = None,
    branch_threshold: float = 0.5,
    time_key: str = "pseudotime",
    # ── Biology-informed velocity bias ────────────────────────────────────
    groupby: Optional[str] = None,
    progenitor_cluster: Optional[str] = None,
    terminal_clusters: Optional[dict] = None,
    bias_strength: float = 0.0,
    # ── Training ──────────────────────────────────────────────────────────
    n_archetypes: int = 5,
    n_epochs: int = 5000,
    key_prefix: str = "scjdo",
    verbose: bool = True,
    **fit_drift_kwargs,
) -> dict:
    """
    Run ``fit_drift`` separately on each branch of a multi-lineage trajectory.

    Avoids the branch-mixing problem where cells from different lineages are
    averaged together in the same pseudotime window. Each branch gets its own
    Jacobian tensor, archetype decomposition, and instability gene ranking.

    Parameters
    ----------
    branch_key  : Key in ``adata.obsm`` (branch probability matrix) or
                  ``adata.obs`` (branch label column).
    branch_names: Names for each branch. If None, reads from
                  ``adata.uns['palantir_branch_names']`` or auto-generates.
    branch_threshold : Minimum branch probability to include a cell (obsm mode).
    time_key    : Pseudotime column (applies to all branches).

    Biology-informed velocity bias
    ------------------------------
    These parameters bias the drift field toward **known biological directions**
    (progenitor → terminal) rather than relying solely on the pseudotime gradient.

    groupby             : ``adata.obs`` column containing cluster labels
                          (e.g. ``'paul15_clusters'``).
    progenitor_cluster  : Cluster label of the root/progenitor population
                          (e.g. ``'7MEP'``).
    terminal_clusters   : Dict mapping branch name → terminal cluster label.
                          Must match the branch names in ``branch_names``.
                          Example: ``{'Ery': '3Ery', 'Neu': '15Mo'}``
    bias_strength       : Strength of the target pull relative to the pseudotime
                          gradient.
                          ``0.0`` — pure pseudotime gradient (default, no bias).
                          ``1.0`` — equal weight (recommended starting point).
                          ``2–5`` — strong bias toward terminal; use when the
                                    branch point is ambiguous.

    Examples
    --------
    **Standard (no bias):**

    >>> models = sjd.tl.fit_drift_branches(adata, branch_key='branch_masks')

    **With known terminal states (Paul15):**

    >>> models = sjd.tl.fit_drift_branches(
    ...     adata,
    ...     branch_key         = 'branch_masks',
    ...     groupby            = 'paul15_clusters',
    ...     progenitor_cluster = '7MEP',
    ...     terminal_clusters  = {'Ery': '3Ery', 'Neu': '15Mo'},
    ...     bias_strength      = 1.5,
    ... )

    **With Palantir branch masks (marrow):**

    >>> models = sjd.tl.fit_drift_branches(
    ...     adata,
    ...     branch_key         = 'branch_masks',
    ...     groupby            = 'cell_type',
    ...     progenitor_cluster = 'HSC',
    ...     terminal_clusters  = {'Ery': 'late_Ery', 'DC': 'DC', 'Mono': 'Mono'},
    ...     bias_strength      = 1.0,
    ... )
    """
    # ── Pre-compute progenitor and terminal centroids for velocity bias ────
    rep = fit_drift_kwargs.get("rep", adata.uns.get("scjdo_prep", {}).get("rep", "X_pca"))
    if rep not in adata.obsm:
        rep = "X_pca"
    X_full = adata.obsm[rep].astype(np.float32)

    progenitor_centroid = None
    if progenitor_cluster is not None and groupby is not None:
        prog_mask = (adata.obs[groupby] == progenitor_cluster).values
        if prog_mask.sum() > 0:
            progenitor_centroid = X_full[prog_mask].mean(0)
            if verbose:
                print(f"[bias] Progenitor '{progenitor_cluster}': "
                      f"{prog_mask.sum()} cells, centroid computed in {rep}")
    # ── Determine branch assignments ───────────────────────────────────────
    if branch_key in adata.obsm:
        bp_raw = adata.obsm[branch_key]

        # Auto-detect branch names: Palantir stores them as column names when
        # save_as_df=True, or in adata.uns['palantir_branch_names']
        if hasattr(bp_raw, "columns"):          # DataFrame (Palantir save_as_df=True)
            auto_names = list(bp_raw.columns)
            bp = bp_raw.values
        else:
            auto_names = (adata.uns.get("palantir_branch_names")
                          or adata.uns.get("branch_names")
                          or [f"branch_{i}" for i in range(bp_raw.shape[1])])
            bp = np.asarray(bp_raw)

        names = branch_names or auto_names

        if bp.dtype == bool or bp.dtype == np.bool_:
            # Palantir boolean masks from select_branch_cells:
            # Each column = all cells on that branch path (including shared progenitors).
            # Biologically correct: progenitors appear in ALL branches because they
            # precede every fate. This preserves the full trajectory per lineage.
            cell_branch = {name: np.where(bp[:, k])[0]
                           for k, name in enumerate(names)}
        else:
            # Float probability matrix — use max-assignment with threshold
            bp = bp.astype(np.float32)
            assignments = np.argmax(bp, axis=1)
            max_probs   = bp.max(axis=1)
            cell_branch = {name: np.where((assignments == k) & (max_probs >= branch_threshold))[0]
                           for k, name in enumerate(names)}

    elif branch_key in adata.obs.columns:
        labels = adata.obs[branch_key].astype(str).values
        names  = branch_names or sorted(set(labels))
        cell_branch = {name: np.where(labels == name)[0] for name in names}

    else:
        raise KeyError(
            f"branch_key='{branch_key}' not found in adata.obsm or adata.obs. "
            f"Run sjd.pp.prepare_trajectory(pseudotime_method='palantir') first, "
            f"or pass a column name from adata.obs."
        )

    models = {}
    for branch_name, idx in cell_branch.items():
        if len(idx) < 50:
            warnings.warn(
                f"Branch '{branch_name}' has only {len(idx)} cells "
                f"(threshold={branch_threshold}). Skipping.",
                UserWarning, stacklevel=2,
            )
            continue

        if verbose:
            print(f"\n{'='*60}")
            print(f"Branch: {branch_name}  ({len(idx)} cells)")
            print(f"{'='*60}")

        # Subset adata to branch cells
        adata_branch = adata[idx].copy()
        key_added    = f"{key_prefix}_{branch_name}"

        # ── Inject biased velocity prior if terminal clusters provided ────
        # For each branch we know the terminal state → bias velocity toward it.
        # This constrains the drift field to follow known biology rather than
        # relying purely on the pseudotime gradient at the noisy branch point.
        if (bias_strength > 0 and terminal_clusters is not None
                and branch_name in terminal_clusters
                and groupby is not None):

            term_cluster = terminal_clusters[branch_name]
            term_mask    = (adata.obs[groupby] == term_cluster).values
            if term_mask.sum() > 0:
                term_centroid = X_full[term_mask].mean(0)

                # Compute biased velocity for branch cells only
                X_branch  = adata_branch.obsm[rep].astype(np.float32)
                pt_branch = adata_branch.obs[time_key].values.astype(np.float32)
                # Progenitor centroid computed in the full dataset
                prog_cen  = progenitor_centroid  # may be None (still works)

                V_biased = _biased_velocity(
                    X_branch, pt_branch,
                    terminal_centroid   = term_centroid,
                    progenitor_centroid = prog_cen,
                    bias_strength       = bias_strength,
                    k                   = 15,
                )
                # Store as a custom layer so fit_drift uses it
                adata_branch.obsm["_biased_velocity"] = V_biased

                if verbose:
                    print(f"  [bias] Terminal '{term_cluster}': "
                          f"{term_mask.sum()} cells | "
                          f"bias_strength={bias_strength}")
                    if prog_cen is not None:
                        axis  = term_centroid - prog_cen
                        align = float(np.dot(V_biased.mean(0), axis) /
                                      (np.linalg.norm(axis) + 1e-8))
                        print(f"  [bias] Alignment with prog→terminal axis: {align:.3f}")

                # Override the velocity prior inside fit_drift via a hook:
                # We temporarily monkey-patch _pseudotime_velocity for this call.
                import scjdo.tl._drift as _drift_mod
                _orig_vel = _drift_mod._pseudotime_velocity

                def _patched_vel(X, pt, k=15):
                    return V_biased   # ignore X, pt — use pre-computed biased vel

                _drift_mod._pseudotime_velocity = _patched_vel
                try:
                    model = fit_drift(
                        adata_branch,
                        time_key     = time_key,
                        n_archetypes = n_archetypes,
                        n_epochs     = n_epochs,
                        key_added    = "scjdo_branch",
                        verbose      = verbose,
                        **fit_drift_kwargs,
                    )
                finally:
                    _drift_mod._pseudotime_velocity = _orig_vel   # always restore
            else:
                warnings.warn(
                    f"Terminal cluster '{term_cluster}' for branch '{branch_name}' "
                    f"not found in adata.obs['{groupby}']. Running without bias.",
                    UserWarning, stacklevel=2,
                )
                model = fit_drift(
                    adata_branch,
                    time_key     = time_key,
                    n_archetypes = n_archetypes,
                    n_epochs     = n_epochs,
                    key_added    = "scjdo_branch",
                    verbose      = verbose,
                    **fit_drift_kwargs,
                )
        else:
            model = fit_drift(
                adata_branch,
                time_key     = time_key,
                n_archetypes = n_archetypes,
                n_epochs     = n_epochs,
                key_added    = "scjdo_branch",
                verbose      = verbose,
                **fit_drift_kwargs,
            )

        # Store branch results back in the full adata
        adata.uns[key_added] = adata_branch.uns["scjdo_branch"]
        adata.uns[key_added]["branch_name"]  = branch_name
        adata.uns[key_added]["branch_cells"] = idx.tolist()
        adata.uns[key_added]["n_cells"]      = int(len(idx))

        # Store per-branch obsm arrays so plotting functions can use them
        # These are stored inside uns (not adata.obsm) to avoid shape conflicts
        if "X_drift" in adata_branch.obsm:
            adata.uns[key_added]["X_drift"]           = adata_branch.obsm["X_drift"]
        if "X_velocity_pseudo" in adata_branch.obsm:
            adata.uns[key_added]["X_velocity_pseudo"] = adata_branch.obsm["X_velocity_pseudo"]
        models[branch_name] = model

        if verbose:
            r2 = adata.uns[key_added]["r2"]
            print(f"  Done. R²={r2:.3f}  stored in adata.uns['{key_added}']")

    if verbose:
        print(f"\nBranch results stored: "
              f"{[f'{key_prefix}_{n}' for n in models]}")

    return models


# ---------------------------------------------------------------------------
# Post-fit per-branch analysis (instability genes + regulators)
# ---------------------------------------------------------------------------

# Canonical columns returned by `infer_regulators`. Used to build an empty
# DataFrame when a branch yields no qualifying TFs, so downstream code can
# treat every branch uniformly.
_EMPTY_REG_COLS = [
    "regulator", "weighted_score", "mean_instability",
    "enrichment_score", "branch_specificity", "peak_archetype",
    "db_confidence", "n_targets", "enrichment_pval", "top_targets",
]


def branch_drift_analysis(
    adata,
    branch_models,
    *,
    key_prefix: str = "scjdo",
    n_genes: int = 15,
    organism: str = "human",
    min_targets: int = 1,
    n_top_regulators: int = 15,
    plot: bool = True,
    save_dir: Optional[str] = None,
    verbose: bool = True,
):
    """
    Run instability-gene + regulator analysis for each branch produced by
    ``fit_drift_branches`` and reconcile results onto the full ``AnnData``.

    Wraps three pieces of bookkeeping that are repeated across notebooks:

    1. **Subset/roundtrip** — copies the branch result onto a per-branch
       cell subset before calling ``pl.instability_genes`` /
       ``tl.infer_regulators`` (mirrors the pattern from the manuscript
       notebooks).
    2. **uns write-back** — ``infer_regulators`` writes to
       ``adata_subset.uns[reg_key]``. This helper copies that entry back onto
       the full ``adata.uns`` under ``f"{key_prefix}_regulators_{branch}"``
       so the ``pl.regulator_*`` plotters can find it.
    3. **Empty-regulator fallback** — catches the ``ValueError`` raised when
       a branch has no qualifying TFs and substitutes an empty DataFrame with
       the canonical column set so callers can iterate without guards.

    Parameters
    ----------
    branch_models : dict
        Output of ``fit_drift_branches``. Only the keys (branch names) are
        used; the values are ignored.
    key_prefix : str
        Must match the ``key_prefix`` passed to ``fit_drift_branches``.
    n_genes : int
        Forwarded to ``pl.instability_genes``.
    organism, min_targets, n_top_regulators
        Forwarded to ``tl.infer_regulators``.
    plot : bool
        If False, calls ``pl.instability_genes`` only when ``save_dir`` is set
        (so PDFs can still be written). When True (default), the function
        produces one figure per branch via the existing plotter.
    save_dir : str, optional
        If given, writes per-branch artifacts:
        ``instab_{branch}.pdf``,
        ``instability_genes_{branch}.csv``,
        ``regulators_{branch}.csv``.
        Directory is created if missing.

    Returns
    -------
    df_genes : dict[str, pandas.DataFrame]
        Branch name → instability-gene table.
    df_regs  : dict[str, pandas.DataFrame]
        Branch name → regulator table (empty DataFrame for branches with no
        qualifying TFs).

    Examples
    --------
    >>> models = sjd.tl.fit_drift_branches(adata, branch_key='branch_masks', ...)
    >>> df_genes, df_regs = sjd.tl.branch_drift_analysis(
    ...     adata, models, organism='human', save_dir='results/figure3_fa/')
    """
    import os
    import pandas as pd

    from scjdo.pl._drift import instability_genes as _plot_instability_genes
    from ._regulators import infer_regulators

    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)

    df_genes: dict = {}
    df_regs:  dict = {}

    for name in branch_models:
        key = f"{key_prefix}_{name}"
        if key not in adata.uns:
            if verbose:
                print(f"[{name}] '{key}' not in adata.uns — skipping")
            continue

        cell_idx = np.asarray(adata.uns[key]["branch_cells"])
        ad_b = adata[cell_idx].copy()
        # `instability_genes` and `infer_regulators` read from
        # ad_b.uns[key]; the entry won't carry over from the subset alone.
        ad_b.uns[key] = adata.uns[key]

        # ── Instability genes ──────────────────────────────────────────────
        instab_save = (os.path.join(save_dir, f"instab_{name}.pdf")
                       if save_dir is not None else None)
        # The plotter always produces a figure and returns the table; we skip
        # calling it only when both `plot=False` and no save path is set.
        if plot or instab_save is not None:
            table = _plot_instability_genes(
                ad_b, key=key, n_genes=n_genes, save=instab_save,
            )
        else:
            table = _plot_instability_genes(
                ad_b, key=key, n_genes=n_genes, save=None,
            )
        df_genes[name] = table
        if save_dir is not None:
            table.to_csv(
                os.path.join(save_dir, f"instability_genes_{name}.csv"),
                index=False,
            )

        # ── Regulators ─────────────────────────────────────────────────────
        reg_key = f"{key_prefix}_regulators_{name}"
        try:
            df_reg = infer_regulators(
                ad_b, key=key, organism=organism,
                min_targets=min_targets, n_top=n_top_regulators,
                key_added=reg_key, verbose=False,
            )
            if reg_key in ad_b.uns:
                adata.uns[reg_key] = ad_b.uns[reg_key]
        except ValueError as e:
            if verbose:
                print(f"  [{name}] infer_regulators: {e} — empty fallback")
            df_reg = pd.DataFrame(columns=_EMPTY_REG_COLS)
        df_regs[name] = df_reg
        if save_dir is not None:
            df_reg.to_csv(
                os.path.join(save_dir, f"regulators_{name}.csv"),
                index=False,
            )

        if verbose:
            top_g = table["gene"].head(10).tolist() if len(table) else []
            top_r = (df_reg["regulator"].head(10).tolist()
                     if len(df_reg) else ["(none)"])
            print(f"\n{name}:")
            print(f"  Top instab genes: {top_g}")
            print(f"  Top regulators:   {top_r}")

    return df_genes, df_regs


# ---------------------------------------------------------------------------
# Interior-peak detection — guard against boundary artifacts in the
# sensitivity profile. The first and last few pseudotime windows have
# systematically lower effective sample size (the kernel mass leaks off the
# end of the trajectory), and the drift field is least constrained there.
# Peaks pinned at t ≈ 0 or t ≈ 1 are usually edge artifacts, not commitment
# events. This helper enforces an interior mask + optional n_eff floor at
# the peak before reporting it.
# ---------------------------------------------------------------------------

def peak_interior(
    res: dict,
    *,
    interior: tuple = (0.10, 0.90),
    n_eff_min_at_peak: Optional[float] = 50.0,
    field: str = "max_real_eig",
) -> dict:
    """
    Extract the *interior-window* peak of a sensitivity curve, guarded
    against boundary artifacts and low-effective-sample-size edges.

    Parameters
    ----------
    res : dict
        A per-branch result dict from :func:`fit_drift` /
        :func:`fit_drift_branches`, i.e. ``adata.uns['scjdo_{branch}']``.
        Must contain ``t_centers`` and the chosen ``field``. If the result
        was produced with kernel windowing, ``n_eff`` is consumed too.
    interior : (lo, hi)
        Pseudotime range considered "interior". The boundary windows
        outside this range are masked out before the peak is computed.
        Default ``(0.10, 0.90)`` is conservative for 100-window scans and
        the 200-grid kernel mode.
    n_eff_min_at_peak : float or None
        If the result carries an ``n_eff`` array (kernel mode) and this is
        set, windows with ``n_eff < n_eff_min_at_peak`` are also masked out.
        Pass ``None`` to disable.
    field : str
        Which curve to peak on. Default ``'max_real_eig'`` (sensitivity);
        any 1-D array in ``res`` is accepted.

    Returns
    -------
    dict
        Always returns a dict with these keys; ``None`` values signal that
        no defensible interior peak exists:

        - ``t``           : pseudotime at the peak, or ``None``
        - ``value``       : peak value of ``field``, or ``None``
        - ``idx``         : grid index of the peak, or ``None``
        - ``interior``    : the ``(lo, hi)`` mask actually applied
        - ``n_valid``     : how many windows survived the mask
        - ``n_eff_at_peak``: ``n_eff`` at the picked window (kernel only)
        - ``boundary_peak``: True if ``argmax`` over the full curve is
                             outside the interior (i.e. the raw report
                             *would* have been a boundary artifact)

    Examples
    --------
    >>> for branch in branch_models:
    ...     res = ad.uns[f'scjdo_{branch}']
    ...     pk  = sjd.tl.peak_interior(res, interior=(0.10, 0.90))
    ...     if pk['t'] is None:
    ...         print(f'  {branch}: no defensible interior peak '
    ...               f'(raw argmax at t={res["t_centers"][res[field].argmax()]:.2f})')
    ...     else:
    ...         flag = ' [WAS BOUNDARY]' if pk['boundary_peak'] else ''
    ...         print(f'  {branch}: t={pk["t"]:.2f}, '
    ...               f'lambda_max={pk["value"]:+.3f}{flag}')
    """
    t   = np.asarray(res["t_centers"], dtype=float)
    eig = np.asarray(res[field],         dtype=float)
    lo, hi = float(interior[0]), float(interior[1])

    raw_argmax = int(np.argmax(eig)) if eig.size else -1
    raw_t      = float(t[raw_argmax]) if raw_argmax >= 0 else None
    boundary   = (raw_t is not None) and ((raw_t < lo) or (raw_t > hi))

    mask = (t >= lo) & (t <= hi)
    n_eff_arr = res.get("n_eff", None)
    if n_eff_arr is not None and n_eff_min_at_peak is not None:
        n_eff_arr = np.asarray(n_eff_arr, dtype=float)
        if n_eff_arr.shape == eig.shape:
            mask &= (n_eff_arr >= float(n_eff_min_at_peak))

    if not mask.any():
        return dict(t=None, value=None, idx=None,
                    interior=(lo, hi), n_valid=0,
                    n_eff_at_peak=None, boundary_peak=boundary)

    masked_eig = np.where(mask, eig, -np.inf)
    idx = int(np.argmax(masked_eig))
    n_eff_at = (float(n_eff_arr[idx])
                if n_eff_arr is not None and n_eff_arr.shape == eig.shape
                else None)
    return dict(t=float(t[idx]), value=float(eig[idx]), idx=idx,
                interior=(lo, hi), n_valid=int(mask.sum()),
                n_eff_at_peak=n_eff_at, boundary_peak=boundary)
