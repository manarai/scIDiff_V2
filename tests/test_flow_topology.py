"""
Tests for scjdo.atlas.flow_topology — FTLE and basin-based fork detection.

Uses a canonical 2D bistable drift v(x,y) = (x - x^3, -y) so the two
attractors are (+1, 0) and (-1, 0), the separatrix is the y-axis, and
the saddle at the origin is a fate boundary. We verify that:
  1. FTLE peaks on the separatrix (near x = 0), not uniformly.
  2. Basin assignment splits seeds by sign(x).
  3. fork_detector confirms only seeds that BOTH ridge AND sit between
     basins — the geometric definition of a fork.
"""
from __future__ import annotations

import numpy as np
import pytest
import torch

from scjdo.atlas.flow_topology import (
    assign_basins,
    compute_ftle,
    fork_detector,
)


class BistableDrift(torch.nn.Module):
    """v(x, y) = (x - x^3, -alpha * y). Two attractors: (+1, 0), (-1, 0)."""

    def __init__(self, alpha: float = 2.0):
        super().__init__()
        self.alpha = alpha

    def forward(self, x: torch.Tensor, t: torch.Tensor) -> torch.Tensor:
        x0 = x[:, 0]
        x1 = x[:, 1]
        return torch.stack([x0 - x0.pow(3), -self.alpha * x1], dim=1)


@pytest.fixture()
def drift():
    return BistableDrift(alpha=2.0)


@pytest.fixture()
def grid_seeds():
    # Skip x=0 exactly: it is the unstable fixed point and stays fixed under
    # deterministic Euler. Any real analysis nudges off it or samples it
    # stochastically; the tests below exercise the mid-strip (|x| < 0.5) via
    # the offset grid instead.
    xs = torch.linspace(-1.4, 1.4, 15) + 0.05
    ys = torch.linspace(-0.4, 0.4, 5)
    gx, gy = torch.meshgrid(xs, ys, indexing="ij")
    return torch.stack([gx.reshape(-1), gy.reshape(-1)], dim=1)


class TestComputeFTLE:
    def test_ftle_higher_near_separatrix(self, drift, grid_seeds):
        ftle = compute_ftle(drift, grid_seeds, t0=0.0, t1=2.0, steps=200)
        assert ftle.shape == (grid_seeds.shape[0],)

        x = grid_seeds[:, 0].numpy()
        near = np.abs(x) < 0.15
        far = np.abs(x) > 1.0
        assert near.sum() > 0 and far.sum() > 0
        assert ftle[near].mean() > ftle[far].mean(), (
            "FTLE should be higher near the x=0 separatrix than deep in "
            "either basin."
        )

    def test_ftle_finite(self, drift, grid_seeds):
        ftle = compute_ftle(drift, grid_seeds, t0=0.0, t1=1.0, steps=100)
        assert np.all(np.isfinite(ftle))


class TestAssignBasins:
    def test_two_basins_by_sign(self, drift, grid_seeds):
        labels, terminal = assign_basins(
            drift, grid_seeds, t0=0.0, t1=4.0, steps=400, n_attractors=2
        )
        assert labels.shape == (grid_seeds.shape[0],)
        assert terminal.shape == grid_seeds.shape

        # Terminal x should be ~+1 or ~-1; label should encode the sign.
        x_end = terminal[:, 0]
        assert np.all(np.abs(np.abs(x_end) - 1.0) < 0.2), (
            f"Terminal x should be near ±1, got range "
            f"[{x_end.min():.3f}, {x_end.max():.3f}]"
        )
        # Same-sign terminal cells should share a label.
        pos_labels = set(labels[x_end > 0.5].tolist())
        neg_labels = set(labels[x_end < -0.5].tolist())
        assert pos_labels.isdisjoint(neg_labels)


class TestForkDetector:
    def test_forks_lie_near_separatrix(self, drift, grid_seeds):
        result = fork_detector(
            drift, grid_seeds,
            t0=0.0, t1=4.0, steps=400,
            n_attractors=2, ftle_percentile=75.0, k_neighbours=6,
        )
        assert result["is_fork"].shape == (grid_seeds.shape[0],)
        assert result["is_fork"].any(), "No forks detected on bistable system"

        x = grid_seeds[:, 0].numpy()
        # All confirmed forks should be within |x| < 0.5 (the mid-strip).
        fork_x = x[result["is_fork"]]
        assert np.all(np.abs(fork_x) < 0.5), (
            f"Fork seeds should sit near the separatrix; got x values "
            f"{sorted(fork_x.tolist())}"
        )

    def test_deep_basin_seeds_not_flagged(self, drift):
        # Seeds well inside each basin, far from the separatrix.
        seeds = torch.tensor(
            [[+1.2, 0.0], [+1.3, 0.1], [-1.2, 0.0], [-1.3, -0.1]],
            dtype=torch.float32,
        )
        # Add a fake near-separatrix seed so ftle_percentile makes sense.
        seeds = torch.cat([seeds, torch.tensor([[0.02, 0.0]], dtype=torch.float32)], dim=0)
        drift = BistableDrift(alpha=2.0)
        result = fork_detector(
            drift, seeds, t0=0.0, t1=4.0, steps=400,
            n_attractors=2, ftle_percentile=80.0, k_neighbours=2,
        )
        # The four deep-basin seeds must not be flagged.
        assert not result["is_fork"][:4].any(), (
            "Deep-basin seeds were incorrectly flagged as forks"
        )
