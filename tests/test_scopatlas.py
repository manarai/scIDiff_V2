"""
Unit tests for scOpAtlas (Stable Operator Atlas).

The old suite stubbed the drift model via unittest.mock, assigning
mock_forward to the instance's __call__ dunder — which Python ignores,
because dunders resolve on the type, not the instance. Every stubbed
model call actually hit the auto-generated stub and returned a fresh
stub object, so no test ever exercised autograd and the compute_all_metrics
no_grad bug went unnoticed. Every test here uses a real nn.Module.
"""

import unittest

import anndata as ad
import numpy as np
import torch
import torch.nn as nn

from scjdo.atlas.atlas_builder import StableOperatorAtlas
from scjdo.atlas.operator_metrics import OperatorMetrics
from scjdo.atlas.regime_classifier import OperatorRegimeClassifier
from scjdo.atlas import statistics as stats_mod
from scjdo.atlas.ensemble import OperatorEnsemble


class _LinearDrift(nn.Module):
    """f(x, t) = x @ A.T. Jacobian is exactly A everywhere, independent of t.

    With this fixture eigenvalues are known in closed form, so metric
    assertions can be tight (5 dp). t is accepted but ignored.
    """

    def __init__(self, A):
        super().__init__()
        self.register_buffer("A", torch.as_tensor(A, dtype=torch.float32))

    def forward(self, x, t):
        return x @ self.A.T


def _diag_drift():
    # Eigenvalues exactly {0.5, -1.0}: one stable, one unstable direction.
    return _LinearDrift(torch.diag(torch.tensor([0.5, -1.0])))


def _rotation_drift():
    # Pure rotation, eigenvalues exactly {+2i, -2i}. Real parts are 0.
    # Needed here (spectrum verified) and by Task 2.2's plasticity-vs-
    # gauge-freedom regression.
    return _LinearDrift([[0.0, -2.0], [2.0, 0.0]])


def _asymmetric_drift():
    # Non-symmetric, non-diagonal Jacobian — guards against silent
    # J vs J.T transposition errors that symmetric/diagonal fixtures hide.
    return _LinearDrift([[1.0, 2.0], [3.0, 4.0]])


class TestOperatorMetrics(unittest.TestCase):
    """Metric computations against a drift whose Jacobian is known."""

    def setUp(self):
        self.drift_model = _diag_drift()
        # A[0,0] = 0.5, A[1,1] = -1.0
        self.expected_A = self.drift_model.A.clone()

        self.metrics_computer = OperatorMetrics(
            drift_model=self.drift_model, epsilon=0.1, device="cpu"
        )

    def test_compute_jacobian(self):
        x = torch.tensor([[1.0, 2.0]])
        t = torch.tensor([0.5])

        jacobian = self.metrics_computer.compute_jacobian(x, t)

        self.assertEqual(jacobian.shape, (1, 2, 2))
        # f_j(x) = sum_k A[j,k] x[k]  =>  d f_j / d x_i = A[j,i]
        # so J should equal A, not A.T. Diagonal fixture would pass either
        # way — the asymmetric fixture (below) is the real transposition
        # regression.
        np.testing.assert_allclose(
            jacobian[0].numpy(), self.expected_A.numpy(), rtol=1e-5, atol=1e-5
        )

    def test_compute_jacobian_asymmetric(self):
        """Regression against silent J vs J.T flips."""
        drift = _asymmetric_drift()
        metrics = OperatorMetrics(drift_model=drift, epsilon=0.1, device="cpu")

        jacobian = metrics.compute_jacobian(
            torch.tensor([[0.7, -0.3]]), torch.tensor([0.5])
        )

        # J must equal A, not A.T. If the convention were flipped, the off-
        # diagonals would swap: [[1,3],[2,4]] instead of [[1,2],[3,4]].
        np.testing.assert_allclose(
            jacobian[0].numpy(), drift.A.numpy(), rtol=1e-5, atol=1e-5
        )

    def test_compute_eigenvalues(self):
        x = torch.tensor([[1.0, 2.0]])
        t = torch.tensor([0.5])

        eigenvalues = self.metrics_computer.compute_eigenvalues(x, t)

        self.assertEqual(eigenvalues.shape, (1, 2))
        eig_real_sorted = np.sort(eigenvalues.real.numpy()[0])
        np.testing.assert_allclose(
            eig_real_sorted, np.sort([-1.0, 0.5]), rtol=1e-5, atol=1e-5
        )

    def test_max_unstable_eigenvalue(self):
        lambda_max = self.metrics_computer.max_unstable_eigenvalue(
            torch.tensor([[1.0, 2.0]]), torch.tensor([0.5])
        )
        self.assertAlmostEqual(lambda_max.item(), 0.5, places=5)

    def test_stability_depth(self):
        lambda_min = self.metrics_computer.stability_depth(
            torch.tensor([[1.0, 2.0]]), torch.tensor([0.5])
        )
        self.assertAlmostEqual(lambda_min.item(), -1.0, places=5)

    def test_plasticity_index(self):
        # With epsilon=0.1, |Re(eig)| in {0.5, 1.0} — neither is near-neutral.
        plasticity = self.metrics_computer.plasticity_index(
            torch.tensor([[1.0, 2.0]]), torch.tensor([0.5])
        )
        self.assertAlmostEqual(plasticity.item(), 0.0, places=5)

    def test_stable_subspace_dim(self):
        stable_dim = self.metrics_computer.stable_subspace_dim(
            torch.tensor([[1.0, 2.0]]), torch.tensor([0.5])
        )
        self.assertEqual(stable_dim.item(), 1.0)

    def test_plasticity_uses_modulus_not_real_part(self):
        """Task 2.1 acceptance regression.

        Rotation Jacobian [[0,-2],[2,0]] has eigenvalues ±2i:
          - |Re(λ)| = 0, 0  → old (buggy) code:   plasticity = 1.0
          - |λ|    = 2, 2   → new (modulus) code: plasticity = 0.0

        Rotation is precisely the gauge-freedom axis (the divergence-free
        circulation the density cannot see). The old |Re| criterion made the
        least identifiable part of the spectrum drive the plasticity signal;
        the modulus definition matches the README and does not.
        """
        drift = _rotation_drift()
        metrics = OperatorMetrics(drift_model=drift, epsilon=0.1, device="cpu")

        p_single = metrics.plasticity_index(
            torch.tensor([[0.3, 0.4]]), torch.tensor([0.5])
        )
        self.assertAlmostEqual(p_single.item(), 0.0, places=5)

        # Same fix must also apply in the batched compute_all_metrics path
        # (there was a duplicated |Re| check there).
        out = metrics.compute_all_metrics(
            torch.randn(6, 2), torch.rand(6), batch_size=3
        )
        np.testing.assert_allclose(out["plasticity"], np.zeros(6), atol=1e-5)

    def test_symmetric_only_uses_symmetric_part(self):
        """Task 2.1: symmetric_only=True eigendecomposes (J + J.T)/2.

        Rotation J = [[0,-2],[2,0]] has (J + J.T)/2 = 0, so all eigenvalues
        of the symmetric part are 0. Plasticity should be 1.0 (all near-
        neutral), stable_dim should be 0 (no strictly-negative eigenvalues).
        This is the gauge-robust reading: pure rotation carries no volume-
        contracting or -expanding structure.
        """
        drift = _rotation_drift()
        metrics = OperatorMetrics(
            drift_model=drift, epsilon=0.1, device="cpu", symmetric_only=True
        )
        eig = metrics.compute_eigenvalues(
            torch.tensor([[0.3, 0.4]]), torch.tensor([0.5])
        )
        np.testing.assert_allclose(eig.real.numpy()[0], [0.0, 0.0], atol=1e-5)

        p = metrics.plasticity_index(
            torch.tensor([[0.3, 0.4]]), torch.tensor([0.5])
        )
        self.assertAlmostEqual(p.item(), 1.0, places=5)

        s = metrics.stable_subspace_dim(
            torch.tensor([[0.3, 0.4]]), torch.tensor([0.5])
        )
        self.assertEqual(s.item(), 0.0)

    def test_rotation_fixture_spectrum(self):
        """Rotation Jacobian [[0,-2],[2,0]] has eigenvalues ±2i.

        Kept here so Task 2.2 (plasticity gauge-invariance) can rely on the
        fixture already being covered.
        """
        drift = _rotation_drift()
        metrics = OperatorMetrics(drift_model=drift, epsilon=0.1, device="cpu")
        eig = metrics.compute_eigenvalues(
            torch.tensor([[0.3, 0.4]]), torch.tensor([0.5])
        )
        # Real parts are zero; magnitudes are 2.
        np.testing.assert_allclose(eig.real.numpy()[0], [0.0, 0.0], atol=1e-5)
        np.testing.assert_allclose(
            np.sort(eig.abs().numpy()[0]), [2.0, 2.0], atol=1e-5
        )

    def test_compute_jacobian_falls_back_when_model_breaks_vmap(self):
        """Regression: real scjdo DriftField uses a kNN velocity module
        that calls .cpu().numpy() during forward — that fails under vmap
        because vmap-batched tensors don't have accessible storage.
        compute_jacobian must fall back to per-cell jacrev in that case
        and still return a correct Jacobian.
        """
        d = 3

        class KNNLikeDrift(nn.Module):
            """Small MLP whose forward triggers the same 'Cannot access
            data pointer' failure the real DriftField hits under vmap."""

            def __init__(self, d):
                super().__init__()
                self.W = nn.Parameter(torch.randn(d, d) * 0.3)
                self.b = nn.Parameter(torch.zeros(d))

            def forward(self, x, t):
                # A branch that materializes a numpy view. Under vmap
                # this raises "Cannot access data pointer of Tensor
                # that doesn't have storage" — matches DriftField's
                # kNN velocity module.
                _ = x.detach().cpu().numpy()
                return x @ self.W.T + self.b

        model = KNNLikeDrift(d).eval()
        X = torch.randn(5, d)
        T = torch.rand(5)

        J_vec = OperatorMetrics(model).compute_jacobian(X, T)

        # Reference: cell-by-cell torch.autograd.functional.jacobian.
        J_ref = torch.stack([
            torch.autograd.functional.jacobian(
                lambda xi, ti=T[i]: model(xi.unsqueeze(0), ti.unsqueeze(0)).squeeze(0),
                X[i],
            )
            for i in range(len(X))
        ])
        self.assertEqual(J_vec.shape, J_ref.shape)
        self.assertLess((J_vec - J_ref).abs().max().item(), 1e-5)

    def test_compute_jacobian_matches_autograd_functional(self):
        """Task 1.3 regression: vectorized jacrev must equal autograd on a
        3-layer asymmetric MLP at dim=8.

        Independent reference is torch.autograd.functional.jacobian applied
        cell-by-cell — different code path from torch.func.jacrev, so a
        transposition or vmap-axis bug would show up here.
        """
        class _MLP(nn.Module):
            def __init__(self, d=8):
                super().__init__()
                # Asymmetric shapes intentional.
                self.net = nn.Sequential(
                    nn.Linear(d + 1, 13), nn.Tanh(),
                    nn.Linear(13, 17), nn.GELU(),
                    nn.Linear(17, d),
                )

            def forward(self, x, t):
                if t.dim() == 1:
                    t = t[:, None]
                return self.net(torch.cat([x, t], dim=1))

        torch.manual_seed(42)
        model = _MLP(8).eval()
        X = torch.randn(6, 8)
        T = torch.rand(6)

        J_vec = OperatorMetrics(model).compute_jacobian(X, T)

        # Cell-by-cell reference via torch.autograd.functional.jacobian.
        J_ref = torch.stack([
            torch.autograd.functional.jacobian(
                lambda xi, ti=T[i]: model(xi.unsqueeze(0), ti.unsqueeze(0)).squeeze(0),
                X[i],
            )
            for i in range(len(X))
        ])
        self.assertEqual(J_vec.shape, J_ref.shape)
        max_diff = (J_vec - J_ref).abs().max().item()
        self.assertLess(max_diff, 1e-5, f"max|J_vec - J_ref| = {max_diff:.2e}")

    def test_compute_all_metrics_end_to_end(self):
        """Regression for the Task 1.1 fix.

        compute_all_metrics used to wrap its batch loop in torch.no_grad(),
        which broke the compute_jacobian call inside it. This test would
        have caught that; the previous stub-based fixture did not, because
        it never actually intercepted the model call and every test errored
        earlier on zeros_like(<stub>).
        """
        drift = _diag_drift()
        metrics_computer = OperatorMetrics(drift_model=drift, epsilon=0.1, device="cpu")

        x = torch.randn(8, 2)
        t = torch.rand(8)

        # Default: eigenvalues NOT returned (Task 1.5); stability_depth_q
        # ADDED (Task 2.2 — quantile-based dimension-stable stability metric).
        out = metrics_computer.compute_all_metrics(x, t, batch_size=4)
        self.assertEqual(
            set(out.keys()),
            {
                "lambda_max_plus",
                "lambda_min_minus",
                "stability_depth_q",
                "plasticity",
                "stable_dim",
            },
        )
        for key in out:
            self.assertEqual(out[key].shape, (8,))

        # Values are exact because the Jacobian is constant.
        np.testing.assert_allclose(out["lambda_max_plus"], np.full(8, 0.5), atol=1e-5)
        np.testing.assert_allclose(out["lambda_min_minus"], np.full(8, -1.0), atol=1e-5)
        np.testing.assert_allclose(out["stable_dim"], np.full(8, 1.0), atol=1e-5)

        # Opt-in: eigenvalues returned as complex64.
        out2 = metrics_computer.compute_all_metrics(x, t, batch_size=4, store_eigenvalues=True)
        self.assertIn("eigenvalues", out2)
        self.assertEqual(out2["eigenvalues"].shape, (8, 2))
        self.assertEqual(out2["eigenvalues"].dtype, np.complex64)


class TestOperatorRegimeClassifier(unittest.TestCase):
    """Classifier tests do not need a drift model — they take a metrics dict."""

    def setUp(self):
        # Legacy-threshold classifier — this class was written against the
        # old hand-set 0.1 / 0.05 / -1.0 / 0.3 thresholds. Opt in explicitly
        # (a) so the tests still document what those numbers did, and (b) so
        # they don't hit the auto path with 1-cell metrics dicts (degenerate:
        # tau_upper would be that single cell's own lambda_max).
        self.classifier = OperatorRegimeClassifier(
            thresholds='legacy',
            plasticity_threshold=0.3,
        )

    def test_classify_stable(self):
        metrics = {
            "lambda_max_plus": np.array([-0.5]),
            "lambda_min_minus": np.array([-0.5]),
            "stability_depth_q": np.array([-0.5]),
            "plasticity": np.array([0.2]),
            "stable_dim": np.array([5.0]),
        }
        regimes, _ = self.classifier.classify(metrics)
        self.assertEqual(regimes[0], "stable")

    def test_classify_plastic(self):
        metrics = {
            "lambda_max_plus": np.array([0.01]),
            "lambda_min_minus": np.array([-0.3]),
            "stability_depth_q": np.array([-0.3]),
            "plasticity": np.array([0.5]),
            "stable_dim": np.array([3.0]),
        }
        regimes, _ = self.classifier.classify(metrics)
        self.assertEqual(regimes[0], "plastic")

    def test_classify_unstable(self):
        metrics = {
            "lambda_max_plus": np.array([0.5]),
            "lambda_min_minus": np.array([-0.2]),
            "stability_depth_q": np.array([-0.2]),
            "plasticity": np.array([0.3]),
            "stable_dim": np.array([2.0]),
        }
        regimes, _ = self.classifier.classify(metrics)
        self.assertEqual(regimes[0], "unstable")

    def test_classify_deeply_stable(self):
        metrics = {
            "lambda_max_plus": np.array([-0.2]),
            "lambda_min_minus": np.array([-2.0]),
            "stability_depth_q": np.array([-2.0]),
            "plasticity": np.array([0.1]),
            "stable_dim": np.array([8.0]),
        }
        regimes, _ = self.classifier.classify(metrics)
        self.assertEqual(regimes[0], "deeply_stable")

    def test_regimes_are_a_partition(self):
        """Task 2.2 regression.

        Four regimes must be mutually exclusive AND cover every cell. Old
        precedence-chain design left this fragile; new construction from two
        orthogonal axes (λ_max band × stability depth) makes it invariant.
        Also confirmed by the internal assertion inside classify().
        """
        metrics = _synth_metrics_from_spectra(n_cells=200, dim=25, seed=1)
        cls = OperatorRegimeClassifier(thresholds="auto")
        _, masks = cls.classify(metrics)
        stacked = np.stack([masks[r] for r in ("stable", "plastic", "unstable", "deeply_stable")])
        np.testing.assert_array_equal(
            stacked.sum(axis=0), np.ones(200, dtype=np.int64)
        )

    def test_hyperbolic_mode_uses_absolute_zero_threshold(self):
        """New default 'hyperbolic' mode: unstable = (λ_max⁺ > 0) absolute.

        A dataset in which every cell has λ_max⁺ < 0 should get 0%
        unstable — NOT 33% as the tercile mode would force. This is the
        load-bearing property that makes cross-dataset regime comparison
        meaningful. Tercile mode ('rank') pins unstable at exactly 33%
        by construction and is retained only for reproducing prior runs.
        """
        rng = np.random.default_rng(4)
        n, dim = 400, 25
        # All cells fully hyperbolic (all eigenvalues negative). Tight
        # variance so max-of-25 stays negative for every cell.
        eigs = rng.normal(loc=-1.5, scale=0.15, size=(n, dim))
        metrics = {
            "lambda_max_plus": eigs.max(axis=1),
            "lambda_min_minus": eigs.min(axis=1),
            "stability_depth_q": np.quantile(eigs, 0.05, axis=1),
            "plasticity": (np.abs(eigs) < 0.1).mean(axis=1),
            "stable_dim": (eigs < 0).sum(axis=1).astype(float),
        }
        # Sanity — this fixture is genuinely all-stable.
        self.assertTrue((metrics["lambda_max_plus"] < 0).all())

        r_hyp, _ = OperatorRegimeClassifier(thresholds="hyperbolic").classify(metrics)
        r_rank, _ = OperatorRegimeClassifier(thresholds="rank").classify(metrics)

        frac_unstable_hyp = float((r_hyp == "unstable").mean())
        frac_unstable_rank = float((r_rank == "unstable").mean())

        # Hyperbolic mode: no unstable cells (correct — no cell has λ_max > 0).
        self.assertLess(frac_unstable_hyp, 0.02,
                        f"hyperbolic: got {frac_unstable_hyp:.1%} unstable")
        # Rank mode: pins at 33% by construction (tercile).
        self.assertGreater(frac_unstable_rank, 0.30,
                           f"rank: got {frac_unstable_rank:.1%} unstable "
                           f"(should be ~33% by tercile construction)")

    def test_dimension_invariance_of_regime_fractions(self):
        """Task 2.2 ACCEPTANCE.

        Spectral distribution is fixed (N(-0.4, 0.5) per eigenvalue); only
        the number of eigenvalues per cell (dim) changes. Under the old
        legacy classifier `min(Re(λ)) < -1.0` swept nearly all non-unstable
        cells into 'deeply_stable' at dim=100 (a minimum-order statistic
        gets more negative as dim grows). Under the new auto classifier,
        regime fractions must stay approximately dim-invariant.
        """
        fractions_new = {}
        fractions_legacy = {}
        for dim in (10, 25, 50, 100):
            m = _synth_metrics_from_spectra(n_cells=400, dim=dim, seed=dim)
            r_new, _ = OperatorRegimeClassifier(thresholds="auto").classify(m)
            r_leg, _ = OperatorRegimeClassifier(thresholds="legacy").classify(m)
            fractions_new[dim] = {name: float((r_new == name).mean()) for name in
                                  ("stable", "plastic", "unstable", "deeply_stable")}
            fractions_legacy[dim] = {name: float((r_leg == name).mean()) for name in
                                     ("stable", "plastic", "unstable", "deeply_stable")}

        # New (auto) classifier: every regime fraction stays tight across
        # dim. Loose 0.10 tolerance because at small dim the tercile split
        # on λ_max jitters a bit; at large dim it's essentially exact.
        for name in ("stable", "plastic", "unstable", "deeply_stable"):
            vals = [fractions_new[d][name] for d in (10, 25, 50, 100)]
            self.assertLess(
                max(vals) - min(vals), 0.10,
                f"AUTO classifier: {name} fraction not dim-invariant, "
                f"got {dict(zip((10, 25, 50, 100), vals))}",
            )

        # Legacy classifier: at least one regime fraction shifts by ≥ 0.15
        # across dim. On this spectrum the shift shows up in 'unstable'
        # (80.5% at dim=10 → 100% at dim=100 — max(Re(λ)) grows with dim, so
        # everything trips the fixed lambda_max > 0.1 unstable threshold),
        # but the concrete mechanism can differ under different spectra;
        # what matters for the acceptance argument is that the legacy path
        # is NOT dim-invariant while the auto path is.
        legacy_max_spread = max(
            max(fractions_legacy[d][name] for d in (10, 25, 50, 100)) -
            min(fractions_legacy[d][name] for d in (10, 25, 50, 100))
            for name in ("stable", "plastic", "unstable", "deeply_stable")
        )
        self.assertGreater(
            legacy_max_spread, 0.15,
            f"legacy path was expected to shift some regime by ≥ 0.15 across "
            f"dim; got max spread {legacy_max_spread:.3f}. Legacy fractions: "
            f"{fractions_legacy}",
        )

    def test_sweep_thresholds_shape(self):
        """Task 2.2 sub-action 4 helper: returns a grid the caller can plot."""
        m = _synth_metrics_from_spectra(n_cells=200, dim=25, seed=2)
        grid = OperatorRegimeClassifier(thresholds="auto").sweep_thresholds(m, n_grid=5)
        self.assertEqual(grid["tau_grid"].shape, (5,))
        self.assertEqual(grid["deep_grid"].shape, (5,))
        for name in ("stable", "plastic", "unstable", "deeply_stable"):
            self.assertEqual(grid[name].shape, (5, 5))
            # Fractions in [0, 1].
            self.assertTrue((grid[name] >= 0).all() and (grid[name] <= 1).all())
        self.assertIn("tau_upper", grid["reference"])

    def test_get_regime_statistics(self):
        regimes = np.array(["stable", "plastic", "stable", "unstable"])
        metrics = {
            "lambda_max_plus": np.array([-0.5, 0.0, -0.3, 0.5]),
            "lambda_min_minus": np.array([-1.0, -0.5, -0.8, -0.2]),
            "stability_depth_q": np.array([-1.0, -0.5, -0.8, -0.2]),
            "plasticity": np.array([0.2, 0.5, 0.3, 0.4]),
            "stable_dim": np.array([5.0, 3.0, 4.0, 2.0]),
        }
        stats = self.classifier.get_regime_statistics(regimes, metrics)
        self.assertEqual(stats["stable"]["count"], 2)
        self.assertAlmostEqual(stats["stable"]["lambda_max_mean"], -0.4, places=5)
        self.assertEqual(stats["plastic"]["count"], 1)
        self.assertAlmostEqual(stats["plastic"]["plasticity_mean"], 0.5, places=5)
        self.assertEqual(stats["unstable"]["count"], 1)
        self.assertAlmostEqual(stats["unstable"]["lambda_max_mean"], 0.5, places=5)


def _synth_metrics_from_spectra(n_cells, dim, seed, epsilon=0.1, q=0.05):
    """Synthesize the metrics dict directly from a prescribed eigenvalue
    distribution — one draw of `dim` real eigenvalues per cell from a fixed
    Gaussian, regardless of ambient dim.

    Bypasses OperatorMetrics/compute_jacobian entirely: the point of the
    cross-dim invariance test is the *classifier*, and we want to isolate
    it from any drift-model dynamics that would themselves change with dim.
    """
    rng = np.random.default_rng(seed)
    eigs = rng.normal(loc=-0.4, scale=0.5, size=(n_cells, dim))
    return {
        "lambda_max_plus": eigs.max(axis=1),
        "lambda_min_minus": eigs.min(axis=1),
        "stability_depth_q": np.quantile(eigs, q, axis=1),
        "plasticity": (np.abs(eigs) < epsilon).mean(axis=1),
        "stable_dim": (eigs < 0).sum(axis=1).astype(float),
    }


def _make_toy_adata(n_cells=60, n_genes=30, n_pcs=4, seed=0):
    """AnnData with the fields StableOperatorAtlas expects."""
    rng = np.random.default_rng(seed)
    adata = ad.AnnData(rng.standard_normal((n_cells, n_genes)).astype(np.float32))
    adata.obsm["X_pca"] = rng.standard_normal((n_cells, n_pcs)).astype(np.float32)
    adata.obs["pseudotime"] = np.linspace(0, 1, n_cells)
    half = n_cells // 2
    adata.obs["cell_type"] = ["TypeA"] * half + ["TypeB"] * (n_cells - half)
    q = n_cells // 4
    adata.obs["condition"] = (
        ["Control"] * q + ["Treatment"] * q + ["Control"] * q + ["Treatment"] * (n_cells - 3 * q)
    )
    return adata


class TestStableOperatorAtlas(unittest.TestCase):
    """Real drift model, real AnnData, real autograd."""

    def setUp(self):
        self.n_pcs = 4
        self.adata = _make_toy_adata(n_pcs=self.n_pcs)
        # 4-dim linear drift with known spectrum {0.5, -0.2, -0.7, -1.5}.
        A = torch.diag(torch.tensor([0.5, -0.2, -0.7, -1.5]))
        self.drift_model = _LinearDrift(A)

    def test_initialization(self):
        atlas = StableOperatorAtlas(
            adata=self.adata,
            drift_model=self.drift_model,
            use_rep="X_pca",
            pseudotime_key="pseudotime",
        )
        self.assertEqual(atlas.use_rep, "X_pca")
        self.assertEqual(atlas.pseudotime_key, "pseudotime")

    def test_validation_missing_rep(self):
        with self.assertRaises(ValueError):
            StableOperatorAtlas(
                adata=self.adata,
                drift_model=self.drift_model,
                use_rep="X_missing",
                pseudotime_key="pseudotime",
            )

    def test_validation_missing_pseudotime_raises(self):
        """Task 1.4 regression.

        Previously this only warned and then substituted np.linspace(0, 1, n),
        assigning pseudotime by AnnData row order — an entirely meaningless
        atlas with a warning that would be scrolled past. Must now raise.
        """
        adata = self.adata.copy()
        del adata.obs["pseudotime"]
        with self.assertRaisesRegex(ValueError, "pseudotime"):
            StableOperatorAtlas(
                adata=adata,
                drift_model=self.drift_model,
                use_rep="X_pca",
                pseudotime_key="pseudotime",
            )

    def test_regime_masks_populated_on_confidence_path(self):
        """Task 1.5 regression.

        build(compute_confidence=True) used to leave self.regime_masks None
        while self.regimes was populated — any downstream code reading masks
        got AttributeError. Masks must now always be present.
        """
        atlas = StableOperatorAtlas(
            adata=self.adata,
            drift_model=self.drift_model,
            use_rep="X_pca",
            pseudotime_key="pseudotime",
        )
        atlas.build(batch_size=8, compute_confidence=True, quality_check=False)
        self.assertIsNotNone(atlas.regime_masks)
        for regime in ("stable", "plastic", "unstable", "deeply_stable"):
            self.assertIn(regime, atlas.regime_masks)
        # Masks partition the cells: exactly one True per row across regimes.
        stacked = np.stack([atlas.regime_masks[r] for r in atlas.regime_masks])
        np.testing.assert_array_equal(stacked.sum(axis=0), np.ones(len(self.adata)))

    def test_boundary_distance_column_named_honestly(self):
        """Task 2.3 regression.

        The confidence path used to write `adata.obs['regime_confidence']`,
        which read as a calibrated posterior probability but was actually a
        clipped, hand-set-threshold-normalized distance to the decision
        boundary — no probabilistic meaning. Column is renamed to
        `boundary_distance` until Task 3.2 replaces it with ensemble
        concordance.
        """
        atlas = StableOperatorAtlas(
            adata=self.adata,
            drift_model=self.drift_model,
            use_rep="X_pca",
            pseudotime_key="pseudotime",
        )
        atlas.build(batch_size=8, compute_confidence=True, quality_check=False)
        self.assertIn("boundary_distance", atlas.adata.obs.columns)
        self.assertNotIn("regime_confidence", atlas.adata.obs.columns)
        vals = atlas.adata.obs["boundary_distance"].to_numpy()
        self.assertTrue(((vals >= 0.0) & (vals <= 1.0)).all())

    def test_store_eigenvalues_default_off(self):
        """Task 1.5 regression.

        Default builds must NOT dump per-cell eigenvalue spectra into
        adata.uns — at 100k × 50 that's ~40 MB written into every saved .h5ad.
        Opt-in via store_eigenvalues=True; opt-in must produce complex64
        (not complex128).
        """
        atlas = StableOperatorAtlas(
            adata=self.adata.copy(),
            drift_model=self.drift_model,
            use_rep="X_pca",
            pseudotime_key="pseudotime",
        )
        atlas.build(batch_size=8, quality_check=False)
        self.assertNotIn("operator_eigenvalues", atlas.adata.uns)
        self.assertNotIn("eigenvalues", atlas.metrics)

        atlas2 = StableOperatorAtlas(
            adata=self.adata.copy(),
            drift_model=self.drift_model,
            use_rep="X_pca",
            pseudotime_key="pseudotime",
        )
        atlas2.build(batch_size=8, store_eigenvalues=True, quality_check=False)
        self.assertIn("operator_eigenvalues", atlas2.adata.uns)
        eig = atlas2.adata.uns["operator_eigenvalues"]
        self.assertEqual(eig.shape, (len(self.adata), self.n_pcs))
        # complex64, not complex128 — halves the on-disk footprint.
        self.assertEqual(eig.dtype, np.complex64)

    def test_end_to_end_build(self):
        """atlas.build() + validate_nonredundancy on a real nn.Module.

        Previously nothing in the module had been run end-to-end on a real
        model; this is that coverage.
        """
        atlas = StableOperatorAtlas(
            adata=self.adata,
            drift_model=self.drift_model,
            use_rep="X_pca",
            pseudotime_key="pseudotime",
        )

        atlas.build(batch_size=8, quality_check=False)

        # Results are stored on the object.
        self.assertIsNotNone(atlas.metrics)
        self.assertIsNotNone(atlas.regimes)
        self.assertEqual(len(atlas.regimes), len(self.adata))

        # …and mirrored into adata.obs.
        for col in (
            "operator_regime",
            "lambda_max_plus",
            "lambda_min_minus",
            "plasticity",
            "stable_dim",
        ):
            self.assertIn(col, atlas.adata.obs.columns)

        # Constant Jacobian ⇒ constant metrics.
        np.testing.assert_allclose(
            atlas.adata.obs["lambda_max_plus"].to_numpy(),
            np.full(len(self.adata), 0.5),
            atol=1e-4,
        )
        np.testing.assert_allclose(
            atlas.adata.obs["lambda_min_minus"].to_numpy(),
            np.full(len(self.adata), -1.5),
            atol=1e-4,
        )

        # Non-redundancy validator runs and returns the expected keys.
        validation = atlas.validate_nonredundancy(
            celltype_key="cell_type", condition_key="condition"
        )
        self.assertIn("celltype_regime_diversity", validation)
        self.assertIn("condition_comparison", validation)


class TestStatistics(unittest.TestCase):
    """Task 3.1 — statistical inference module."""

    # ---------------- unit tests ----------------

    def test_chi_square_returns_expected_keys(self):
        rng = np.random.default_rng(0)
        n = 200
        regimes = rng.choice(["stable", "unstable", "plastic", "deeply_stable"], size=n)
        conditions = rng.choice(["A", "B"], size=n)
        out = stats_mod.chi_square_regime_vs_condition(regimes, conditions)
        for key in ("chi2", "p", "dof", "cramers_v", "n_cells", "table"):
            self.assertIn(key, out)
        self.assertGreaterEqual(out["cramers_v"], 0.0)
        self.assertLessEqual(out["cramers_v"], 1.0)

    def test_bh_correction_matches_definition(self):
        # From Benjamini & Hochberg (1995), Table 1 sanity: ordered
        # p-values [0.01, 0.02, 0.03, 0.5] with m=4 → q-values increase.
        q = stats_mod.bh_correction([0.01, 0.02, 0.03, 0.5])
        # Monotone non-decreasing under ranking.
        self.assertLess(q[0], q[3])
        self.assertLessEqual(q[0], q[1])
        self.assertLessEqual(q[1], q[2])
        # Never > 1.
        self.assertTrue((q <= 1.0).all())

    def test_bootstrap_ci_shape_and_bounds(self):
        rng = np.random.default_rng(0)
        regimes = rng.choice(["stable", "plastic", "unstable", "deeply_stable"], size=200)
        boot = stats_mod.bootstrap_ci_regime_fractions(regimes, n_boot=200, rng=rng)
        for r, (point, lo, hi) in boot.items():
            self.assertGreaterEqual(lo, 0.0)
            self.assertLessEqual(hi, 1.0)
            self.assertLessEqual(lo, point)
            self.assertLessEqual(point, hi)

    def test_depth_matched_subsample_equalizes_counts(self):
        rng = np.random.default_rng(0)
        # Overlapping-but-shifted depth distributions. If they were nearly
        # disjoint, quantile-bin matching has nothing to match against — a
        # separate documented limitation, not what this test is guarding.
        ct = np.array(["T"] * 400)
        cond = np.array(["A"] * 300 + ["B"] * 100)
        counts = np.concatenate([
            rng.normal(loc=5000, scale=1500, size=300),
            rng.normal(loc=3500, scale=1500, size=100),
        ])
        keep = stats_mod.depth_matched_subsample(ct, cond, counts, n_bins=5, rng=rng)
        # Post-matching, per-condition N is equal within this celltype.
        kept_conds = cond[keep]
        n_A = (kept_conds == "A").sum()
        n_B = (kept_conds == "B").sum()
        self.assertEqual(n_A, n_B)
        # Depth-mean difference must SHRINK (not merely stay flat).
        diff_pre = abs(counts[cond == "A"].mean() - counts[cond == "B"].mean())
        diff_post = abs(counts[keep & (cond == "A")].mean()
                        - counts[keep & (cond == "B")].mean())
        self.assertLess(diff_post, diff_pre * 0.5)

    def test_miller_madow_reduces_bias_at_small_n(self):
        # Two bins, uniform underlying distribution. Plug-in H is
        # downward-biased at small n; MM adds a positive correction that
        # gets it closer to log2(2) = 1.0. The correction is 1/(2n·ln 2),
        # so at n=1000 it's ~0.0007 above the plug-in — nonzero even at
        # "large" n. That is by design; don't assert exact equality to 1.
        counts_small = np.array([3, 3])  # small n
        counts_large = np.array([500, 500])  # large n
        h_small = stats_mod.miller_madow_entropy(counts_small)
        h_large = stats_mod.miller_madow_entropy(counts_large)
        # MM correction adds a positive term so it's ≥ plug-in.
        p = counts_small / counts_small.sum()
        plugin_small = float(-np.sum(p * np.log2(p)))
        self.assertGreater(h_small, plugin_small)  # strict at small n
        # Large-n estimate rides very close to 1.0.
        self.assertAlmostEqual(h_large, 1.0, places=2)
        # And the correction is smaller at large n than at small n
        # (bias correction shrinks with 1/n).
        self.assertLess(h_large - 1.0, h_small - plugin_small)

    # ---------------- integration & acceptance ----------------

    def test_permutation_detects_real_effect(self):
        """Positive control: labels strongly correlated with regime → tiny p."""
        rng = np.random.default_rng(1)
        n = 100
        # Condition A is mostly stable, condition B is mostly unstable.
        regimes = np.array(["stable"] * n + ["unstable"] * n)
        conditions = np.array(["A"] * n + ["B"] * n)
        # Add some noise so table isn't degenerate.
        flip = rng.random(2 * n) < 0.1
        conditions = np.where(flip, np.where(conditions == "A", "B", "A"), conditions)
        res = stats_mod.permutation_p_value(regimes, conditions, n_perm=500, rng=rng)
        self.assertLess(res["p_perm"], 0.01)
        self.assertGreater(res["cramers_v"], 0.5)

    def test_ACCEPTANCE_negative_control(self):
        """Task 3.1 ACCEPTANCE — must pass before anything else in Phase 3.

        Split a single homogeneous cohort into two arbitrary halves,
        label them 'A' and 'B', run the full pipeline. Median p-value
        across seeds must stay well above 0.05. If a random split
        produces significant non-redundancy, the test is measuring
        drift-field noise, not biology, and the paper's headline claim
        has no support.

        We run 15 independent seeded splits so a single unlucky draw
        cannot fail the test (under H0 ~5% of draws will nominally
        reject at p=0.05; the MEDIAN is the stable summary).
        """
        rng_master = np.random.default_rng(0)
        n = 300  # per celltype
        # One celltype, homogeneous regime distribution.
        celltype = np.array(["T"] * n)
        # Regime and count distributions the same for every cell — the
        # split is the only source of variation.
        regimes = rng_master.choice(
            ["stable", "plastic", "unstable", "deeply_stable"],
            size=n, p=[0.5, 0.2, 0.15, 0.15],
        )
        counts = rng_master.normal(loc=5000, scale=800, size=n).clip(min=100)

        p_values = []
        matched_medians = []
        for seed in range(15):
            rng = np.random.default_rng(seed + 100)
            # Random 50/50 split.
            condition = rng.choice(["A", "B"], size=n)
            res = stats_mod.test_nonredundancy(
                regimes, celltype, condition, counts=counts,
                n_perm=200, n_matched=15, n_perm_matched=100, n_boot=200,
                rng=rng,
            )
            p_values.append(res.per_celltype["T"]["p_perm"])
            matched_medians.append(res.matched["T"]["p_perm_median"])

        # Median observed p must be well away from 0.05 under H0.
        self.assertGreater(
            float(np.median(p_values)), 0.15,
            f"acceptance failed: median p_perm across 15 random splits was "
            f"{np.median(p_values):.3f} ≤ 0.15 — the test is picking up "
            f"noise on homogeneous data. p_values: {p_values}",
        )
        # And matched p-values should behave the same way.
        self.assertGreater(
            float(np.median(matched_medians)), 0.15,
            f"acceptance failed: median matched p across 15 splits was "
            f"{np.median(matched_medians):.3f} ≤ 0.15",
        )

    def test_test_nonredundancy_smoke(self):
        """End-to-end: real effect present in one celltype, none in another.

        Confirms the result dict has the promised shape and that BH
        q-values are computed across cell types.
        """
        rng = np.random.default_rng(2)
        # 2 celltypes, 2 conditions, one celltype has an effect and one doesn't.
        n = 150
        celltype = np.array(["Effect"] * n + ["Null"] * n)
        # Effect celltype: A→mostly stable, B→mostly unstable.
        cond_eff = rng.choice(["A", "B"], size=n)
        r_eff = np.where(cond_eff == "A", "stable", "unstable")
        # 10% flip to keep it messy.
        flip = rng.random(n) < 0.1
        r_eff = np.where(flip, np.where(r_eff == "stable", "unstable", "stable"), r_eff)
        # Null celltype: uniform regimes, random conditions.
        cond_null = rng.choice(["A", "B"], size=n)
        r_null = rng.choice(["stable", "plastic", "unstable", "deeply_stable"], size=n)

        regimes = np.concatenate([r_eff, r_null])
        conditions = np.concatenate([cond_eff, cond_null])
        counts = rng.normal(loc=5000, scale=800, size=2 * n).clip(min=100)

        res = stats_mod.test_nonredundancy(
            regimes, celltype, conditions, counts=counts,
            n_perm=300, n_matched=20, n_perm_matched=100, n_boot=200, rng=rng,
        )

        # Shape.
        self.assertIn("Effect", res.per_celltype)
        self.assertIn("Null", res.per_celltype)
        self.assertIn("Effect", res.q_values)
        self.assertIn("Null", res.q_values)
        self.assertTrue(res.depth_matching_applied)

        # Effect celltype: significant, large effect size.
        self.assertLess(res.per_celltype["Effect"]["p_perm"], 0.05)
        self.assertGreater(res.per_celltype["Effect"]["cramers_v"], 0.4)
        # Null celltype: not significant.
        self.assertGreater(res.per_celltype["Null"]["p_perm"], 0.05)


class _RandomInitMLPDrift(nn.Module):
    """Small MLP drift model whose weights are drawn from a torch RNG.

    Used to construct ensembles where members differ only in random
    initialization (weakest form of ensemble diversity — real ensembles
    should also vary velocity-prior and other structural axes, but this
    is enough to exercise the concordance calculation).
    """

    def __init__(self, dim, seed):
        super().__init__()
        gen = torch.Generator().manual_seed(seed)
        self.net = nn.Sequential(nn.Linear(dim + 1, 16), nn.Tanh(), nn.Linear(16, dim))
        for p in self.parameters():
            p.data = torch.randn(p.shape, generator=gen) * 0.3

    def forward(self, x, t):
        if t.dim() == 1:
            t = t[:, None]
        return self.net(torch.cat([x, t], dim=1))


class TestOperatorEnsemble(unittest.TestCase):
    """Task 3.2 — ensemble across drift-field members."""

    def setUp(self):
        self.n_pcs = 4
        self.adata = _make_toy_adata(n_pcs=self.n_pcs, n_cells=80)

    def test_requires_at_least_two_members(self):
        with self.assertRaises(ValueError):
            OperatorEnsemble([_diag_drift()])

    def test_identical_members_give_perfect_concordance(self):
        """Task 3.2 sanity: if every member is the SAME model, concordance
        must be exactly 1.0 for every cell — no ensemble noise possible."""
        m = _RandomInitMLPDrift(self.n_pcs, seed=42)
        ens = OperatorEnsemble([m, m, m])
        ens.build(self.adata, verbose=False)
        self.assertTrue((ens.concordance == 1.0).all())
        self.assertEqual(ens.discordant_fraction(0.5), 0.0)

    def test_different_members_produce_meaningful_concordance(self):
        """Task 3.2 acceptance-flavored test: an ensemble of independently
        initialized MLP drifts must produce a non-trivial concordance
        distribution and write both columns into adata.obs.
        """
        adata = self.adata.copy()
        members = [_RandomInitMLPDrift(self.n_pcs, seed=s) for s in (1, 2, 3, 4, 5)]
        ens = OperatorEnsemble(members)
        ens.build(adata, verbose=False)

        self.assertEqual(ens.n_members, 5)
        # concordance in the [1/N, 1] interval.
        self.assertGreaterEqual(ens.concordance.min(), 1 / 5 - 1e-9)
        self.assertLessEqual(ens.concordance.max(), 1.0 + 1e-9)
        # Non-trivial spread — if this ever collapses to 1.0 across the
        # whole population we're not exercising the ensemble at all.
        self.assertLess(ens.concordance.min(), 1.0)
        # adata.obs is populated.
        self.assertIn("ensemble_regime", adata.obs.columns)
        self.assertIn("regime_concordance", adata.obs.columns)
        # Report the discordant-fraction headline — expected around 10–30%
        # on this synthetic setup. Guarded loosely: this test doesn't
        # need to pin an exact number, only that a discordance signal
        # exists.
        self.assertGreater(ens.discordant_fraction(1.0), 0.0)

    def test_concordant_mask(self):
        members = [_RandomInitMLPDrift(self.n_pcs, seed=s) for s in (1, 2, 3)]
        ens = OperatorEnsemble(members)
        ens.build(self.adata, verbose=False)
        mask = ens.concordant_mask(threshold=0.75)
        self.assertEqual(mask.shape, (len(self.adata),))
        # threshold 0 keeps everything, threshold >1 keeps nothing.
        self.assertTrue(ens.concordant_mask(threshold=0.0).all())
        self.assertFalse(ens.concordant_mask(threshold=1.001).any())

    def test_concordance_key_wires_into_validate_nonredundancy(self):
        """Task 3.2 + 3.1 integration: with regime_concordance in obs,
        validate_nonredundancy runs a second 'concordant-only' inference
        pipeline and stores the results alongside the all-cells version.
        """
        adata = self.adata.copy()
        drift = _diag_drift_4d = _LinearDrift(
            torch.diag(torch.tensor([0.5, -0.2, -0.7, -1.5]))
        )
        atlas = StableOperatorAtlas(
            adata=adata, drift_model=drift, use_rep="X_pca",
            pseudotime_key="pseudotime",
        )
        atlas.build(batch_size=8, thresholds="legacy", quality_check=False)

        # Fake a concordance column such that ~half the cells are "concordant".
        rng = np.random.default_rng(0)
        adata.obs["regime_concordance"] = rng.uniform(0, 1, size=len(adata))

        val = atlas.validate_nonredundancy(
            celltype_key="cell_type",
            condition_key="condition",
            n_perm=50, n_matched=5, n_perm_matched=25, n_boot=50,
            concordance_threshold=0.5,
            verbose=False,
        )

        # Layer 3 must appear alongside the layer-2 statistics.
        self.assertIn("statistics", val)
        self.assertIn("statistics_concordant_only", val)
        c = val["statistics_concordant_only"]
        self.assertEqual(c["concordance_key"], "regime_concordance")
        self.assertEqual(c["concordance_threshold"], 0.5)
        # Roughly half the cells retained (uniform draw > 0.5).
        self.assertGreater(c["n_cells_kept"], 0)
        self.assertLess(c["n_cells_kept"], c["n_cells_total"])


class TestDriftQualityGate(unittest.TestCase):
    """Task 3.3 — drift-quality gate."""

    def _trajectory_adata(self, n_cells=200, n_pcs=4, seed=0):
        """AnnData whose cells follow a linear trajectory in latent space.

        Cell x_i at pseudotime t_i is drawn from N(t_i · μ, σ²I) with a
        fixed direction μ. A drift model that predicts f(x) = μ will move
        cells in the right direction and improve transport SWD over the
        "stay put" baseline — legitimate trajectory data for the gate.
        The quality gate resorts internally, so we do not need to shuffle
        row order to exercise that path.
        """
        rng = np.random.default_rng(seed)
        pseudotime = np.linspace(0, 1, n_cells)
        mu = np.zeros(n_pcs); mu[0] = 3.0  # trajectory along first PC
        X = np.stack([
            rng.normal(loc=t * mu, scale=0.3, size=n_pcs) for t in pseudotime
        ]).astype(np.float32)
        adata = ad.AnnData(rng.standard_normal((n_cells, 20)).astype(np.float32))
        adata.obsm["X_pca"] = X
        adata.obs["pseudotime"] = pseudotime
        adata.obs["cell_type"] = ["TypeA"] * (n_cells // 2) + \
                                 ["TypeB"] * (n_cells - n_cells // 2)
        return adata, mu

    def test_trained_drift_passes_gate(self):
        """A drift that actually predicts the trajectory beats the
        no-op baseline and passes the transport check."""
        adata, mu = self._trajectory_adata(seed=1)
        # "Trained" drift: constant f(x) = mu. Perfect for this trajectory.
        class ConstantDrift(nn.Module):
            def __init__(self, mu):
                super().__init__()
                self.register_buffer("mu", torch.as_tensor(mu, dtype=torch.float32))
            def forward(self, x, t):
                return self.mu.expand_as(x)
        model = ConstantDrift(mu)

        from scjdo.atlas.quality import check_drift_quality
        report = check_drift_quality(
            model, adata,
            use_rep="X_pca", pseudotime_key="pseudotime",
        )
        self.assertTrue(report["passed"], msg=report["summary"])
        te = report["criteria"]["transport_error"]
        self.assertTrue(te["ran"] and te["passed"])
        # Transport ratio should be substantially < 1 (drift beats baseline).
        self.assertLess(te["ratio"], 0.9)

    def test_ACCEPTANCE_untrained_drift_fails_gate(self):
        """Task 3.3 ACCEPTANCE. A deliberately un-trained drift (random
        MLP output, no training) must FAIL the transport-error check —
        pushing cells forward with random directions can't beat the
        stay-put baseline. If this test ever passes without touching the
        gate, the gate isn't actually gating anything.
        """
        adata, _ = self._trajectory_adata(seed=2)
        # "Under-trained" model: random-init MLP, no training.
        # Small init so predictions are near-noise around zero.
        class UntrainedDrift(nn.Module):
            def __init__(self, d, seed):
                super().__init__()
                gen = torch.Generator().manual_seed(seed)
                self.net = nn.Sequential(nn.Linear(d + 1, 16), nn.Tanh(),
                                          nn.Linear(16, d))
                for p in self.parameters():
                    p.data = torch.randn(p.shape, generator=gen) * 0.5
            def forward(self, x, t):
                if t.dim() == 1: t = t[:, None]
                return self.net(torch.cat([x, t], dim=1))
        model = UntrainedDrift(4, seed=42)

        from scjdo.atlas.quality import check_drift_quality
        report = check_drift_quality(
            model, adata,
            use_rep="X_pca", pseudotime_key="pseudotime",
        )
        self.assertFalse(report["passed"], msg=report["summary"])
        self.assertFalse(report["criteria"]["transport_error"]["passed"])

    def test_gate_is_enforced_in_build_by_default(self):
        """StableOperatorAtlas.build() raises DriftQualityError when the
        gate fails, unless force=True is set."""
        from scjdo.atlas.quality import DriftQualityError
        adata, _ = self._trajectory_adata(seed=3)
        # Untrained drift → gate fails.
        class UntrainedDrift(nn.Module):
            def __init__(self, d, seed):
                super().__init__()
                gen = torch.Generator().manual_seed(seed)
                self.net = nn.Sequential(nn.Linear(d + 1, 16), nn.Tanh(),
                                          nn.Linear(16, d))
                for p in self.parameters():
                    p.data = torch.randn(p.shape, generator=gen) * 0.5
            def forward(self, x, t):
                if t.dim() == 1: t = t[:, None]
                return self.net(torch.cat([x, t], dim=1))

        atlas = StableOperatorAtlas(
            adata=adata, drift_model=UntrainedDrift(4, seed=99),
            use_rep="X_pca", pseudotime_key="pseudotime",
        )
        with self.assertRaises(DriftQualityError):
            atlas.build(batch_size=8)
        # Report was still written to uns before the raise.
        self.assertIn("drift_quality", atlas.adata.uns)
        self.assertFalse(atlas.adata.uns["drift_quality"]["passed"])

    def test_force_bypasses_gate_failure(self):
        """force=True proceeds past a failing gate (with a loud banner)."""
        adata, _ = self._trajectory_adata(seed=4)
        class UntrainedDrift(nn.Module):
            def __init__(self, d, seed):
                super().__init__()
                gen = torch.Generator().manual_seed(seed)
                self.net = nn.Sequential(nn.Linear(d + 1, 16), nn.Tanh(),
                                          nn.Linear(16, d))
                for p in self.parameters():
                    p.data = torch.randn(p.shape, generator=gen) * 0.5
            def forward(self, x, t):
                if t.dim() == 1: t = t[:, None]
                return self.net(torch.cat([x, t], dim=1))

        atlas = StableOperatorAtlas(
            adata=adata, drift_model=UntrainedDrift(4, seed=1),
            use_rep="X_pca", pseudotime_key="pseudotime",
        )
        # Should proceed and complete despite failing quality.
        atlas.build(batch_size=8, force=True)
        self.assertIsNotNone(atlas.regimes)
        # And the failure is recorded in adata.uns so downstream code
        # (or a reviewer) can see the gate was overridden.
        self.assertFalse(atlas.adata.uns["drift_quality"]["passed"])

    def test_velocity_check_skipped_when_no_velocity(self):
        """No velocity in adata → velocity criterion 'ran=False', skipped."""
        adata, mu = self._trajectory_adata(seed=5)
        class ConstantDrift(nn.Module):
            def __init__(self, mu):
                super().__init__()
                self.register_buffer("mu", torch.as_tensor(mu, dtype=torch.float32))
            def forward(self, x, t):
                return self.mu.expand_as(x)
        from scjdo.atlas.quality import check_drift_quality
        rep = check_drift_quality(
            ConstantDrift(mu), adata,
            use_rep="X_pca", pseudotime_key="pseudotime",
        )
        self.assertFalse(rep["criteria"]["velocity_agreement"]["ran"])
        self.assertFalse(rep["criteria"]["ensemble_concordance"]["ran"])
        # But transport ran → overall passed (>=1 criterion ran, all passed).
        self.assertTrue(rep["passed"])

    def test_ensemble_criterion_reads_concordance_column(self):
        """When adata.obs['regime_concordance'] exists, the ensemble
        criterion runs and its pass/fail depends on the discordant
        fraction threshold."""
        adata, mu = self._trajectory_adata(seed=6)
        # Case 1: concordance mostly high → passes.
        rng = np.random.default_rng(0)
        adata.obs["regime_concordance"] = rng.uniform(0.7, 1.0, size=len(adata))
        class ConstantDrift(nn.Module):
            def __init__(self, mu):
                super().__init__()
                self.register_buffer("mu", torch.as_tensor(mu, dtype=torch.float32))
            def forward(self, x, t):
                return self.mu.expand_as(x)
        from scjdo.atlas.quality import check_drift_quality
        rep = check_drift_quality(
            ConstantDrift(mu), adata,
            use_rep="X_pca", pseudotime_key="pseudotime",
        )
        self.assertTrue(rep["criteria"]["ensemble_concordance"]["ran"])
        self.assertTrue(rep["criteria"]["ensemble_concordance"]["passed"])

        # Case 2: concordance mostly low → ensemble criterion fails.
        adata.obs["regime_concordance"] = rng.uniform(0.0, 0.3, size=len(adata))
        rep2 = check_drift_quality(
            ConstantDrift(mu), adata,
            use_rep="X_pca", pseudotime_key="pseudotime",
        )
        self.assertFalse(rep2["criteria"]["ensemble_concordance"]["passed"])
        self.assertFalse(rep2["passed"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
