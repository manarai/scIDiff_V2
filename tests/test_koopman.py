"""
Tests for scjdo.archetypes.koopman.

Covers:
  - Return-shape parity with jacobian_modes: (K, D, D), (T, K), float.
  - Non-negativity of activations (magnitudes of complex projections).
  - Recovery of a known growth rate on a synthetic autonomous linear
    system (J(t) = A constant) — the top continuous-time eigenvalue
    should match Re(eig(A)) up to discretisation error.
  - Recovery of oscillation on a rotating linear system.
  - Local mode returns a per-t eigenvalue array of the right shape.
  - decompose_archetypes dispatcher works for both backends.
"""
from __future__ import annotations

import numpy as np
import torch

from scjdo.archetypes.koopman import koopman_modes
from scjdo.archetypes import decompose_archetypes


def _linear_autonomous_tensor(A: np.ndarray, T: int = 60, dt: float = 0.02,
                              noise: float = 0.0, seed: int = 0) -> tuple[torch.Tensor, np.ndarray]:
    """
    Build a Jacobian tensor for a system whose linearisation is a
    constant A but whose observed operator sequence is A perturbed by
    small noise (so EDMD sees a real trajectory instead of a degenerate
    one). Returns (J_tensor, grid).
    """
    rng = np.random.default_rng(seed)
    D = A.shape[0]
    grid = np.arange(T, dtype=np.float64) * dt
    Js = np.tile(A[None, :, :], (T, 1, 1)).astype(np.float32)
    if noise > 0:
        # A slow drift plus jitter so consecutive Jacobians differ but
        # stay close to A — this mimics the smooth pseudotime evolution
        # that scJDO produces.
        drift = noise * np.cumsum(rng.standard_normal((T, D, D)).astype(np.float32), axis=0) / np.sqrt(T)
        Js = Js + drift
    return torch.from_numpy(Js), grid


class TestKoopmanShapes:
    def test_local_mode_shapes(self):
        A = np.array([[0.3, 1.0], [-1.0, 0.3]], dtype=np.float32)   # spiral, Re=0.3
        J, grid = _linear_autonomous_tensor(A, T=40, dt=0.05, noise=0.02, seed=1)
        patterns, activations, err = koopman_modes(
            J, rank=2, mode="local", window_half=5, grid=grid,
        )
        assert patterns.shape == (2, 2, 2)
        assert activations.shape == (40, 2)
        assert (activations.numpy() >= 0).all()
        assert isinstance(err, float) and np.isfinite(err)

    def test_global_mode_shapes(self):
        A = np.array([[0.5, 0.0], [0.0, -0.2]], dtype=np.float32)
        J, grid = _linear_autonomous_tensor(A, T=30, dt=0.05, noise=0.01, seed=2)
        patterns, activations, err = koopman_modes(
            J, rank=2, mode="global", grid=grid,
        )
        assert patterns.shape[1:] == (2, 2)
        assert activations.shape[0] == 30
        assert (activations.numpy() >= 0).all()

    def test_diagnostics_returned(self):
        A = np.array([[0.5, 0.0], [0.0, -0.2]], dtype=np.float32)
        J, grid = _linear_autonomous_tensor(A, T=30, dt=0.05, noise=0.01, seed=3)
        patterns, activations, err, diag = koopman_modes(
            J, rank=2, mode="local", window_half=4, grid=grid,
            return_diagnostics=True,
        )
        required = {"eigenvalues", "growth_rates", "freqs", "V_r",
                    "singular_values", "K_seq", "mode", "reduced_rank",
                    "delta_tau", "nyquist_ang_freq", "nyquist_cycle_freq",
                    "geometry", "whiten"}
        assert required <= set(diag.keys())
        assert diag["eigenvalues"].shape == (30, patterns.shape[0])
        assert diag["growth_rates"].shape == (30, patterns.shape[0])
        assert diag["freqs"].shape == (30, patterns.shape[0])
        assert diag["mode"] == "local"


class TestKoopmanBranchAmbiguityCompanions:
    """Nyquist limit and delta_tau must travel with the numbers."""

    def test_nyquist_matches_delta_tau(self):
        A = np.array([[0.5, 0.0], [0.0, -0.2]], dtype=np.float32)
        dt = 0.05
        J, grid = _linear_autonomous_tensor(A, T=30, dt=dt, noise=0.01, seed=4)
        _, _, _, diag = koopman_modes(
            J, rank=2, mode="local", window_half=4, grid=grid,
            return_diagnostics=True,
        )
        # Nyquist in the same units as `freqs` (rad / unit-τ) is π / Δτ.
        interior = slice(5, 25)   # skip boundary windows with clipped support
        ang = np.asarray(diag["nyquist_ang_freq"])[interior]
        dtau = np.asarray(diag["delta_tau"])[interior]
        assert np.allclose(ang, np.pi / dtau, rtol=1e-6)
        cyc = np.asarray(diag["nyquist_cycle_freq"])[interior]
        assert np.allclose(cyc, 1.0 / (2.0 * dtau), rtol=1e-6)
        # And the interior windows should reflect the sampling step.
        assert np.allclose(dtau, dt, atol=1e-6)

    def test_global_mode_nyquist_is_scalar_broadcast(self):
        A = np.array([[0.5, 0.0], [0.0, -0.2]], dtype=np.float32)
        J, grid = _linear_autonomous_tensor(A, T=30, dt=0.05, noise=0.01, seed=5)
        _, _, _, diag = koopman_modes(
            J, rank=2, mode="global", grid=grid, return_diagnostics=True,
        )
        ang = np.asarray(diag["nyquist_ang_freq"])
        # Every window shares the same Δτ under global mode.
        assert np.allclose(ang, ang[0], atol=1e-10)


class TestKoopmanSpectralRecovery:
    """
    The vectorised trajectory of J(t)=A across pseudotime lives on a
    single vector, so DMD/EDMD there recovers eigenvalue 1 trivially and
    is not the right invariant to test. Instead, test that the *reduced
    trajectory* recovers the correct dominant growth mode when we drive
    the tensor with an exponentially growing amplitude a(t) e^{μ t} A —
    i.e. the trajectory itself follows the mode's evolution.
    """

    def _driven_tensor(self, mu: complex, T: int = 80, dt: float = 0.05,
                       D: int = 3, seed: int = 0) -> tuple[torch.Tensor, np.ndarray]:
        """
        Build J(t) whose vec-trajectory has a controllable spectrum.

        - Real mu (imag == 0): scalar exp(mu t) amplitude on one direction.
          The reduced trajectory is 1-D and DMD recovers exp(mu dt) exactly.
        - Purely imaginary mu (real == 0): two-dimensional rotation
          J(t) = cos(omega t) A1 + sin(omega t) A2 with A1 ⟂ A2 in the
          Frobenius sense. DMD in the resulting 2-D reduced trajectory
          recovers a conjugate pair of eigenvalues at angle omega*dt.
        """
        rng = np.random.default_rng(seed)
        grid = np.arange(T, dtype=np.float64) * dt
        if mu.imag == 0:
            A0 = rng.standard_normal((D, D)).astype(np.float32)
            A0 /= np.linalg.norm(A0)
            amp = np.exp(mu.real * grid).astype(np.float32)
            Js = amp[:, None, None] * A0[None, :, :]
        else:
            A1 = rng.standard_normal((D, D)).astype(np.float32)
            A2 = rng.standard_normal((D, D)).astype(np.float32)
            # Frobenius-orthogonalise A2 against A1 so vec(A1) ⟂ vec(A2).
            inner = float((A1 * A2).sum() / max((A1 * A1).sum(), 1e-8))
            A2 = A2 - inner * A1
            A1 /= np.linalg.norm(A1)
            A2 /= np.linalg.norm(A2) + 1e-8
            envelope = np.exp(mu.real * grid).astype(np.float32)
            c = (envelope * np.cos(mu.imag * grid)).astype(np.float32)
            s = (envelope * np.sin(mu.imag * grid)).astype(np.float32)
            Js = c[:, None, None] * A1[None, :, :] + s[:, None, None] * A2[None, :, :]
        # Tiny jitter so the trajectory has more than one direction after
        # SVD reduction (EDMD needs r >= 1 with variance) even in the
        # real-mu case where the ideal trajectory is 1-D.
        Js = Js + 1e-3 * rng.standard_normal(Js.shape).astype(np.float32)
        return torch.from_numpy(Js.astype(np.float32)), grid

    def test_growth_rate_recovered(self):
        mu_true = 1.5
        J, grid = self._driven_tensor(complex(mu_true, 0.0), T=100, dt=0.02, seed=10)
        _, _, _, diag = koopman_modes(
            J, rank=3, mode="global", grid=grid, return_diagnostics=True,
        )
        # ref_eigenvalues is |λ|-sorted, so the dominant mode is index 0.
        # For a driven-growth trajectory that mode's continuous-time rate
        # should be ~mu_true.
        dt = float(np.median(np.diff(grid)))
        cts = np.log(diag["ref_eigenvalues"]) / dt
        top_re = float(cts[0].real)
        assert abs(top_re - mu_true) < 0.15, f"expected ~{mu_true}, got {top_re}"

    def test_oscillation_frequency_recovered(self):
        omega_true = 5.0
        J, grid = self._driven_tensor(complex(0.0, omega_true), T=200, dt=0.02, seed=11)
        _, _, _, diag = koopman_modes(
            J, rank=3, mode="global", grid=grid, return_diagnostics=True,
        )
        dt = float(np.median(np.diff(grid)))
        cts = np.log(diag["ref_eigenvalues"]) / dt
        # Dominant conjugate pair: entries 0 and 1 of the |λ|-sorted list.
        # (Small-|λ| noise modes trail behind and are ignored here.)
        dominant_omegas = np.abs(cts[:2].imag)
        assert np.all(np.abs(dominant_omegas - omega_true) < 0.5), (
            f"expected ~{omega_true}, got {dominant_omegas}"
        )


class TestKoopmanGeometry:
    """Branch-free descriptors — Henrici / reactivity / transient gain."""

    def _run(self, J, grid, mode="global"):
        _, _, _, diag = koopman_modes(
            J, rank=2, mode=mode, grid=grid, return_diagnostics=True,
        )
        return diag["geometry"]

    def test_normal_operator_has_zero_henrici_and_unit_reactivity(self):
        # Build a Jacobian trajectory whose vec-sequence, after SVD
        # reduction, lands in an orthogonal 2-D subspace so the fitted
        # K is a normal (rotation-like) operator. Two Frobenius-
        # orthogonal matrices driven by (cos, sin) achieve this.
        rng = np.random.default_rng(42)
        D, T, dt, omega = 3, 200, 0.02, 3.0
        A1 = rng.standard_normal((D, D)).astype(np.float32)
        A2 = rng.standard_normal((D, D)).astype(np.float32)
        A2 = A2 - float((A1 * A2).sum() / max((A1 * A1).sum(), 1e-8)) * A1
        A1 /= np.linalg.norm(A1); A2 /= np.linalg.norm(A2) + 1e-8
        grid = np.arange(T) * dt
        c = np.cos(omega * grid).astype(np.float32)
        s = np.sin(omega * grid).astype(np.float32)
        Js = c[:, None, None] * A1[None, :, :] + s[:, None, None] * A2[None, :, :]
        # NO noise — we want the ideal rotation so the geometry
        # descriptors hit their theoretical limits.
        g = self._run(torch.from_numpy(Js.astype(np.float32)), grid, mode="global")
        # Global mode → scalars broadcast along T; sample the first entry.
        assert g["henrici"][0] < 1e-3, f"expected ≈0, got {g['henrici'][0]}"
        assert abs(g["reactivity"][0] - 1.0) < 1e-3, f"expected ≈1, got {g['reactivity'][0]}"

    def test_reactivity_lower_bound(self):
        # Reactivity σ_max(K) / ρ(K) is ≥ 1 for every K by construction.
        rng = np.random.default_rng(7)
        D, T = 4, 40
        Js = rng.standard_normal((T, D, D)).astype(np.float32)
        grid = np.arange(T) * 0.05
        g = self._run(torch.from_numpy(Js), grid, mode="local")
        # Every window ≥ 1 within fp tolerance.
        assert np.all(np.asarray(g["reactivity"]) >= 1.0 - 1e-6)

    def test_transient_gain_nonnegative_and_finite(self):
        rng = np.random.default_rng(8)
        D, T = 4, 30
        Js = rng.standard_normal((T, D, D)).astype(np.float32)
        grid = np.arange(T) * 0.05
        g = self._run(torch.from_numpy(Js), grid, mode="local")
        arr = np.asarray(g["transient_gain"])
        assert np.isfinite(arr).all()
        assert (arr >= 0).all()


class TestKoopmanWhitenMetric:
    """
    Whitened metric changes the eigenvectors (so geometry descriptors
    change), but leaves the eigenvalues invariant — the two operators
    are similar.
    """

    def test_eigenvalues_invariant_under_whitening_at_ridge_zero(self):
        # Under exact LS (ridge=0) the two operators are similar (K and
        # Σ⁻¹ K Σ) so their eigenvalue spectra coincide. With ridge > 0
        # the regulariser acts differently in the two metrics and the
        # eigenvalues shift — that behaviour is documented in the
        # koopman_modes docstring and in KOOPMAN.md §2.
        A = np.array([[0.6, -0.4], [0.4, 0.6]], dtype=np.float32)
        J, grid = _linear_autonomous_tensor(A, T=60, dt=0.03, noise=0.02, seed=9)
        _, _, _, d0 = koopman_modes(J, rank=2, mode="global", grid=grid,
                                    ridge=0.0, whiten=False,
                                    return_diagnostics=True)
        _, _, _, d1 = koopman_modes(J, rank=2, mode="global", grid=grid,
                                    ridge=0.0, whiten=True,
                                    return_diagnostics=True)
        lam0 = np.sort(np.abs(np.asarray(d0["ref_eigenvalues"])))
        lam1 = np.sort(np.abs(np.asarray(d1["ref_eigenvalues"])))
        assert np.allclose(lam0, lam1, atol=1e-4), (
            f"|λ| should match under exact LS: {lam0} vs {lam1}"
        )
        assert d1["whiten"] is True and d0["whiten"] is False


class TestDispatcher:
    def test_snmf_backend_via_dispatcher(self):
        torch.manual_seed(0)
        J = torch.randn(12, 4, 4)
        p, a, err = decompose_archetypes(J, method="snmf", rank=3, n_restarts=1)
        assert p.shape == (3, 4, 4)
        assert a.shape == (12, 3)

    def test_koopman_backend_via_dispatcher(self):
        torch.manual_seed(0)
        J = torch.randn(20, 4, 4)
        out = decompose_archetypes(J, method="koopman", rank=3)
        # Dispatcher sets return_diagnostics=True by default for koopman.
        assert len(out) == 4
        p, a, err, diag = out
        assert p.shape[0] <= 3 and p.shape[1:] == (4, 4)
        assert a.shape[0] == 20
        assert "growth_rates" in diag

    def test_unknown_method_raises(self):
        import pytest
        J = torch.randn(10, 3, 3)
        with pytest.raises(ValueError, match="Unknown method"):
            decompose_archetypes(J, method="nope")


class TestFitDriftIntegration:
    """
    Smoke test: fit_drift with archetype_method='koopman' returns a
    populated adata.uns entry containing a 'koopman' sub-dict with the
    spectral diagnostics. Uses a very small synthetic AnnData so the
    test runs fast.
    """

    def _tiny_adata(self, N: int = 120, D: int = 6, G: int = 20, seed: int = 0):
        import anndata as ad
        rng = np.random.default_rng(seed)
        X = rng.standard_normal((N, G)).astype(np.float32)
        pca = rng.standard_normal((N, D)).astype(np.float32)
        pt = np.linspace(0.0, 1.0, N).astype(np.float32)
        adata = ad.AnnData(X=X)
        adata.obsm["X_pca"] = pca
        adata.obs["pseudotime"] = pt
        adata.varm["PCs"] = rng.standard_normal((G, D)).astype(np.float32)
        return adata

    def test_koopman_integration(self):
        from scjdo.tl import fit_drift
        adata = self._tiny_adata()
        fit_drift(
            adata,
            time_key="pseudotime",
            rep="X_pca",
            n_epochs=25,                 # tiny — just to smoke-test the path
            n_archetypes=3,
            windowing="kernel",
            grid_size=40,
            bandwidth=0.05,
            n_boot=2,
            archetype_method="koopman",
            koopman_kwargs=dict(window_half=4, ridge=1e-3, mode="local"),
            key_added="scjdo_koop",
            verbose=False,
        )
        res = adata.uns["scjdo_koop"]
        assert "patterns" in res and "activations" in res
        assert "koopman" in res, "Koopman diagnostics not stored"
        koop = res["koopman"]
        assert "eigenvalues" in koop and "growth_rates" in koop and "freqs" in koop
        assert res["params"]["archetype_method"] == "koopman"
        # Activations must be non-negative for downstream instability code.
        assert (res["activations"] >= 0).all()
