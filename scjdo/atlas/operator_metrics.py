"""
Operator Metrics Module

Computes operator-derived state metrics from Jacobian eigenvalues.
These metrics define the dynamical properties of cellular states.
"""

import torch
import numpy as np
from torch.func import jacrev, vmap
from typing import Dict, Optional, Tuple


def _compute_jacobian_autograd(drift_model, x, t):
    """Row-by-row Jacobian via classic torch.autograd.grad.

    Used as the fallback when torch.func.vmap/jacrev can't trace through
    the drift model (e.g. because it calls .numpy() internally, as
    scJDO's DriftField does for its kNN velocity prior). Slower than
    the vmap path but robust to any autograd-compatible forward.
    """
    batch_size, dim = x.shape
    jacobians = []
    for i in range(batch_size):
        xi = x[i:i + 1].clone().requires_grad_(True)
        ti = t[i:i + 1]
        drift = drift_model(xi, ti)
        rows = []
        for j in range(dim):
            g = torch.zeros_like(drift)
            g[0, j] = 1.0
            grad = torch.autograd.grad(
                drift, xi, grad_outputs=g,
                retain_graph=True, create_graph=False,
            )[0]
            rows.append(grad)
        jacobians.append(torch.cat(rows, dim=0))
    return torch.stack(jacobians, dim=0)


class OperatorMetrics:
    """
    Compute operator-derived state metrics from Jacobian eigenvalues.

    Given a drift field f(x,t), this class computes the Jacobian J(x,t) = ∂f/∂x
    and extracts four metrics from its eigenvalue spectrum:

    1. λ_max⁺: Max real part of the spectrum
    2. λ_min⁻: Stability depth (min real part)
    3. P: Plasticity index — fraction of eigenvalues with |λ| < epsilon
       (complex modulus, NOT |Re(λ)|)
    4. S: Stable subspace dimension — # of eigenvalues with Re(λ) < 0

    Gauge freedom (IMPORTANT). Snapshot single-cell data determines the drift
    field only up to a divergence-free (rotational) circulation. Consequences:

    - trace(J) = ∇·f is gauge-invariant.
    - # of eigenvalues with Re(λ) < 0 (stability depth, stable_dim) is
      approximately gauge-invariant; sign of Re(λ) is topologically robust.
    - INDIVIDUAL complex eigenvalue pairs are NOT gauge-invariant — rotation
      is precisely the ambiguity, so anything counting purely imaginary /
      near-imaginary eigenvalues from a single drift model has no data support.
    - λ_max⁺ as a *bifurcation detector* is gauge-sensitive at the crossing
      itself. Its sign is robust; its magnitude is not.

    This is why plasticity is defined on |λ| (complex modulus), not |Re(λ)|:
    the old |Re(λ)| criterion counted pure rotations (|λ| large, Re(λ) = 0)
    as maximally plastic, which is the least identifiable part of the spectrum.
    Set symmetric_only=True to eigendecompose (J + J.T)/2 instead — the
    symmetric part governs local volume change and is far more robust to the
    antisymmetric ambiguity; its eigenvalues are real, so the modulus-vs-real-
    part distinction disappears entirely.

    Args:
        drift_model: Trained DriftField model
        epsilon: Threshold for near-neutral modes, in modulus (default: 0.1)
        device: Device for computation (default: "cpu")
        symmetric_only: If True, eigendecompose (J + J.T)/2 instead of J.
            Trades a full spectral picture for gauge robustness. Recommended
            when the goal is stability depth / plasticity rather than
            bifurcation location.
    """

    def __init__(
        self,
        drift_model,
        epsilon: float = 0.1,
        device: str = "cpu",
        symmetric_only: bool = False,
        stability_quantile: float = 0.05,
    ):
        self.drift_model = drift_model
        self.epsilon = epsilon
        self.device = device
        self.symmetric_only = symmetric_only
        # Per-cell quantile of Re(λ) used as a dimension-stable proxy for
        # stability depth. The extremum min(Re(λ)) is a minimum-order
        # statistic whose expected value drifts with dim (in 50 PCs almost
        # every trained model has *some* eigenvalue below -1). A quantile
        # doesn't. Reported in metrics as 'stability_depth_q'.
        self.stability_quantile = stability_quantile
        self.drift_model.to(device)
        self.drift_model.eval()
    
    def compute_jacobian(
        self,
        x: torch.Tensor,
        t: torch.Tensor
    ) -> torch.Tensor:
        """
        Compute Jacobian matrix ∂f/∂x at (x,t).

        Args:
            x: State tensor (batch_size, dim)
            t: Time tensor (batch_size,)

        Returns:
            Jacobian tensor (batch_size, dim, dim), with J[b, j, i] = ∂f_j/∂x_i
            at cell b — the standard mathematical convention that torch.func.jacrev
            returns.
        """
        x = x.to(self.device)
        t = t.to(self.device)

        drift_model = self.drift_model

        def _f_single(xi, ti):
            # xi: (dim,), ti: scalar tensor. Batch to (1, ...) for the model,
            # unbatch its output back to (dim,) so jacrev sees a vector-in,
            # vector-out function.
            return drift_model(xi.unsqueeze(0), ti.unsqueeze(0)).squeeze(0)

        # Fast path: vmap-batched jacrev — ~80× faster on plain nn.Module.
        # Falls back to classic torch.autograd.grad if the drift model
        # contains ops that require tensors with real storage. In
        # particular, scJDO's DriftField uses a kNN velocity lookup that
        # calls .detach().cpu().numpy() during forward — under torch.func's
        # batched tracing (both vmap AND jacrev), tensors are BatchedTensor
        # wrappers with no accessible storage, so that op raises
        # "Cannot access data pointer …". Falling all the way back to
        # eager autograd is the only path that survives.
        try:
            return vmap(jacrev(_f_single))(x, t)
        except RuntimeError as e:
            if "data pointer" not in str(e) and "vmap" not in str(e).lower():
                raise
            return _compute_jacobian_autograd(drift_model, x, t)
    
    def compute_eigenvalues(
        self,
        x: torch.Tensor,
        t: torch.Tensor
    ) -> torch.Tensor:
        """
        Compute eigenvalues of J(x,t) — or of (J + J.T)/2 if symmetric_only.

        Args:
            x: State tensor (batch_size, dim)
            t: Time tensor (batch_size,)

        Returns:
            Eigenvalues (batch_size, dim), complex-valued dtype (complex64 for
            symmetric_only, upcast for consistency with the general path).
        """
        jacobians = self.compute_jacobian(x, t)

        if self.symmetric_only:
            J_sym = 0.5 * (jacobians + jacobians.transpose(-1, -2))
            # Real symmetric ⇒ real eigenvalues via eigvalsh. Cast to complex
            # so downstream metric code (which reads .real / .abs()) is
            # uniform across the two modes.
            eigenvalues = torch.linalg.eigvalsh(J_sym).to(torch.complex64)
        else:
            eigenvalues = torch.linalg.eigvals(jacobians)

        # compute_jacobian requires autograd, but eigenvalues themselves are
        # leaf data for all downstream metric reductions — detach so the graph
        # is released as soon as the batch is done.
        return eigenvalues.detach()
    
    def max_unstable_eigenvalue(
        self,
        x: torch.Tensor,
        t: torch.Tensor
    ) -> torch.Tensor:
        """
        Compute λ_max⁺ = max(Re(λᵢ)).

        The SIGN of λ_max⁺ (stable vs unstable) is topologically robust and
        approximately gauge-invariant. Its MAGNITUDE at bifurcation, and the
        exact spatial location of the crossing, are NOT gauge-invariant —
        interpret with the ensemble-concordance quantifier (Task 3.2) rather
        than as a precise eigenvalue.

        Args:
            x: State tensor (batch_size, dim)
            t: Time tensor (batch_size,)

        Returns:
            Max unstable eigenvalue (batch_size,)
        """
        eigenvalues = self.compute_eigenvalues(x, t)
        real_parts = eigenvalues.real
        return real_parts.max(dim=1)[0]
    
    def stability_depth(
        self,
        x: torch.Tensor,
        t: torch.Tensor
    ) -> torch.Tensor:
        """
        Compute λ_min⁻ = min(Re(λᵢ)).
        
        Measures how strongly deviations are damped.
        More negative values indicate deeper commitment.
        
        Args:
            x: State tensor (batch_size, dim)
            t: Time tensor (batch_size,)
            
        Returns:
            Stability depth (batch_size,)
        """
        eigenvalues = self.compute_eigenvalues(x, t)
        real_parts = eigenvalues.real
        return real_parts.min(dim=1)[0]
    
    def plasticity_index(
        self,
        x: torch.Tensor,
        t: torch.Tensor
    ) -> torch.Tensor:
        """
        Compute P = #{|λᵢ| < ε} / d — complex modulus, matching the README.

        Old code counted #{|Re(λᵢ)| < ε}, which classified pure rotational
        modes (Re(λ) = 0, |λ| large) as maximally plastic. Rotation is
        precisely the gauge-freedom axis, so the old definition made the
        LEAST identifiable part of the spectrum drive the "plasticity"
        signal. Fixed to use |λ| here and in compute_all_metrics.

        Args:
            x: State tensor (batch_size, dim)
            t: Time tensor (batch_size,)

        Returns:
            Plasticity index (batch_size,) in [0, 1]
        """
        eigenvalues = self.compute_eigenvalues(x, t)
        dim = eigenvalues.shape[1]

        near_neutral = (eigenvalues.abs() < self.epsilon).sum(dim=1).float()
        return near_neutral / dim
    
    def stable_subspace_dim(
        self,
        x: torch.Tensor,
        t: torch.Tensor
    ) -> torch.Tensor:
        """
        Compute S = #{Re(λᵢ) < 0}.
        
        Number of stable directions (buffering capacity).
        Higher values indicate more robust states.
        
        Args:
            x: State tensor (batch_size, dim)
            t: Time tensor (batch_size,)
            
        Returns:
            Stable subspace dimension (batch_size,)
        """
        eigenvalues = self.compute_eigenvalues(x, t)
        real_parts = eigenvalues.real
        return (real_parts < 0).sum(dim=1).float()
    
    def compute_all_metrics(
        self,
        x: torch.Tensor,
        t: torch.Tensor,
        batch_size: int = 32,
        store_eigenvalues: bool = False,
    ) -> Dict[str, np.ndarray]:
        """
        Compute all operator metrics for a set of cells.

        Processes data in batches to avoid memory issues.

        Args:
            x: State tensor (n_cells, dim)
            t: Time tensor (n_cells,)
            batch_size: Batch size for processing
            store_eigenvalues: If True, return the full per-cell eigenvalue
                spectra (as complex64) under key 'eigenvalues'. Off by default
                because the spectrum is ~80 MB at 100k cells × 50 dims and
                lands in adata.uns → every saved .h5ad.

        Returns:
            Dictionary with keys:
                - 'lambda_max_plus': Max Re(λ) per cell
                - 'lambda_min_minus': Min Re(λ) per cell (extremum, kept for
                    diagnostics — prefer 'stability_depth_q' for classification)
                - 'stability_depth_q': Per-cell quantile of Re(λ) at
                    self.stability_quantile (default 5th percentile). Dimension-
                    stable proxy for stability depth; use in place of
                    lambda_min_minus when comparing across different dim.
                - 'plasticity': Fraction of |λ| < epsilon per cell
                - 'stable_dim': # of Re(λ) < 0 per cell
                - 'eigenvalues': Full eigenvalue spectra (only if
                    store_eigenvalues=True)
        """
        n_cells = x.shape[0]

        lambda_max_plus = []
        lambda_min_minus = []
        stability_depth_q = []
        plasticity = []
        stable_dim = []
        all_eigenvalues = [] if store_eigenvalues else None

        # compute_jacobian requires autograd, so we cannot wrap this loop in
        # torch.no_grad(). compute_eigenvalues detaches its output, so the
        # per-batch graph is released as soon as eigenvalues are returned.
        for i in range(0, n_cells, batch_size):
            batch_x = x[i:i+batch_size]
            batch_t = t[i:i+batch_size]

            eigenvalues = self.compute_eigenvalues(batch_x, batch_t)

            if store_eigenvalues:
                all_eigenvalues.append(eigenvalues.to(torch.complex64).cpu())

            real_parts = eigenvalues.real
            moduli = eigenvalues.abs()

            lambda_max_plus.append(real_parts.max(dim=1)[0].cpu())
            lambda_min_minus.append(real_parts.min(dim=1)[0].cpu())
            stability_depth_q.append(
                torch.quantile(real_parts, self.stability_quantile, dim=1).cpu()
            )

            dim = real_parts.shape[1]
            # Plasticity: complex modulus, not |Re|. See plasticity_index.
            near_neutral = (moduli < self.epsilon).sum(dim=1).float()
            plasticity.append((near_neutral / dim).cpu())

            stable_dim.append((real_parts < 0).sum(dim=1).float().cpu())

            # Drop the batch spectrum eagerly if we're not keeping it.
            del eigenvalues, real_parts, moduli

        metrics = {
            'lambda_max_plus': torch.cat(lambda_max_plus).numpy(),
            'lambda_min_minus': torch.cat(lambda_min_minus).numpy(),
            'stability_depth_q': torch.cat(stability_depth_q).numpy(),
            'plasticity': torch.cat(plasticity).numpy(),
            'stable_dim': torch.cat(stable_dim).numpy(),
        }
        if store_eigenvalues:
            metrics['eigenvalues'] = torch.cat(all_eigenvalues).numpy()

        return metrics

    def compute_metrics_at_pseudotime(
        self,
        x: torch.Tensor,
        pseudotime: np.ndarray,
        batch_size: int = 32,
        store_eigenvalues: bool = False,
    ) -> Dict[str, np.ndarray]:
        """
        Compute operator metrics using pseudotime as the time coordinate.

        Args:
            x: State tensor (n_cells, dim)
            pseudotime: Pseudotime values (n_cells,)
            batch_size: Batch size for processing
            store_eigenvalues: See compute_all_metrics.

        Returns:
            Dictionary of operator metrics
        """
        # Normalize pseudotime to [0, 1]
        ptime_min = pseudotime.min()
        ptime_max = pseudotime.max()
        t_normalized = (pseudotime - ptime_min) / (ptime_max - ptime_min + 1e-8)

        t = torch.tensor(t_normalized, dtype=torch.float32)

        return self.compute_all_metrics(
            x, t, batch_size=batch_size, store_eigenvalues=store_eigenvalues
        )
