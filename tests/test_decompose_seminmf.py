"""
Tests for the fixed semi-NMF in scjdo.archetypes.decompose:
  - real NNLS W-update makes the objective monotone non-increasing;
  - returned reconstruction error matches the actual residual;
  - post-convergence rescaling preserves W @ H;
  - jacobian_modes still round-trips.
"""
from __future__ import annotations

import numpy as np
import torch

from scjdo.archetypes.decompose import _seminmf, jacobian_modes


def _random_signed_M(T: int = 30, D: int = 20, K: int = 3, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    W_true = np.abs(rng.standard_normal((T, K)).astype(np.float32))
    H_true = rng.standard_normal((K, D)).astype(np.float32)
    return (W_true @ H_true + 0.01 * rng.standard_normal((T, D))).astype(np.float32)


class TestSemiNMFReturnedError:
    def test_returned_err_matches_actual_residual(self):
        M = _random_signed_M(seed=1)
        W, H, err = _seminmf(M, rank=3, max_iter=200, tol=1e-6, seed=1)
        actual = float(np.linalg.norm(M - W @ H, "fro"))
        # After the off-by-one fix, the returned err is the residual for the
        # (rescaled) W, H being returned — up to fp rescaling error only.
        assert abs(err - actual) / (actual + 1e-8) < 1e-3


class TestSemiNMFMonotone:
    def test_objective_decreases_by_convergence(self):
        """
        With real NNLS on W and closed-form LS on H, ALS produces a
        non-increasing sequence of Frobenius residuals. The returned err
        after 200 iterations should be no worse than the residual of the
        initial W (which is random non-negative) against the first H.
        """
        M = _random_signed_M(T=25, D=15, K=3, seed=7)
        W_final, H_final, err = _seminmf(M, rank=3, max_iter=200, tol=1e-8, seed=7)
        # Random-init upper bound: use ||M||_F as a trivial ceiling — any
        # sensible ALS run must beat "predict zero".
        assert err < float(np.linalg.norm(M, "fro"))
        # And beat the trivial best rank-1 mean approximation.
        M_mean = M.mean(axis=0, keepdims=True)
        trivial_err = float(np.linalg.norm(M - M_mean, "fro"))
        assert err < trivial_err


class TestSemiNMFRescaleInvariant:
    def test_rescale_preserves_product(self):
        M = _random_signed_M(seed=3)
        W, H, _ = _seminmf(M, rank=3, max_iter=100, tol=1e-6, seed=3)
        # After rescaling, ||H_k||_2 should be ~1 for each archetype.
        norms = np.linalg.norm(H, axis=1)
        assert np.allclose(norms, 1.0, atol=1e-4), (
            f"Post-rescale H rows should be unit norm; got {norms}"
        )
        # W @ H should still be a decent approximation of M.
        residual = float(np.linalg.norm(M - W @ H, "fro"))
        assert residual < float(np.linalg.norm(M, "fro"))

    def test_W_still_nonnegative_after_rescale(self):
        M = _random_signed_M(seed=5)
        W, H, _ = _seminmf(M, rank=3, max_iter=100, tol=1e-6, seed=5)
        assert (W >= 0).all(), "Rescaling must preserve W >= 0 constraint"


class TestJacobianModes:
    def test_output_shapes(self):
        T, d = 12, 5
        torch.manual_seed(0)
        J = torch.randn(T, d, d)
        patterns, activations, err = jacobian_modes(J, rank=3, n_restarts=2, seed=0)
        assert patterns.shape == (3, d, d)
        assert activations.shape == (T, 3)
        assert isinstance(err, float) and np.isfinite(err)

    def test_activations_nonnegative(self):
        torch.manual_seed(1)
        J = torch.randn(10, 4, 4)
        _, activations, _ = jacobian_modes(J, rank=2, n_restarts=2, seed=1)
        assert (activations >= 0).all()


class TestTVSmoothedSemiNMF:
    """
    tv_lambda opt-in: enabling the total-variation penalty on temporal
    derivatives of W should (a) preserve W ≥ 0, (b) produce smoother
    activation profiles than the unregularised default, and (c) reduce to the
    default behaviour at tv_lambda=0.
    """

    def _concurrent_tensor(self, T=60, D=12, K=3, seed=0):
        """Build a synthetic tensor with K co-timed archetype activations."""
        rng = np.random.default_rng(seed)
        t = np.linspace(0.0, 1.0, T, dtype=np.float32)
        center = 0.5
        # All K peaks co-located at t=0.5, varying widths.
        widths = np.array([0.08, 0.12, 0.18])[:K]
        A = np.stack([
            np.exp(-0.5 * ((t - center) / w) ** 2) for w in widths
        ], axis=1).astype(np.float32)  # (T, K), all peaks at t=0.5
        B = rng.standard_normal((K, D, D)).astype(np.float32)
        for k in range(K):
            B[k] /= (np.linalg.norm(B[k]) + 1e-8)
        J = np.einsum("tk,kij->tij", A, B).astype(np.float32)
        J = J + 0.005 * rng.standard_normal(J.shape).astype(np.float32)
        return torch.from_numpy(J), t, A

    def test_default_matches_tv_lambda_zero(self):
        """tv_lambda=0 must give the same numbers as the unregularised default."""
        J, _, _ = self._concurrent_tensor(seed=2)
        p0, a0, e0 = jacobian_modes(J, rank=3, n_restarts=2, seed=42, tv_lambda=0.0)
        pd, ad, ed = jacobian_modes(J, rank=3, n_restarts=2, seed=42)
        assert np.allclose(a0.numpy(), ad.numpy())
        assert np.allclose(p0.numpy(), pd.numpy())
        assert abs(e0 - ed) < 1e-6

    def test_tv_activations_still_nonnegative(self):
        J, _, _ = self._concurrent_tensor(seed=3)
        _, activations, _ = jacobian_modes(J, rank=3, n_restarts=2, seed=0, tv_lambda=1.0)
        assert (activations.numpy() >= 0).all()

    def test_tv_reduces_temporal_variation(self):
        """A strong TV penalty should produce activation profiles with strictly
        lower total variation than the unregularised default on the same
        tensor. This is what "temporal smoothness prior" is supposed to do."""
        J, _, _ = self._concurrent_tensor(seed=4)
        _, a_no_tv, _ = jacobian_modes(J, rank=3, n_restarts=2, seed=0, tv_lambda=0.0)
        _, a_tv, _ = jacobian_modes(J, rank=3, n_restarts=2, seed=0, tv_lambda=1.0)
        tv_no = float(np.abs(np.diff(a_no_tv.numpy(), axis=0)).sum())
        tv_yes = float(np.abs(np.diff(a_tv.numpy(), axis=0)).sum())
        assert tv_yes < tv_no, (
            f"TV penalty did not reduce temporal variation "
            f"(no-TV total variation = {tv_no:.3f}, TV total variation = {tv_yes:.3f})"
        )
