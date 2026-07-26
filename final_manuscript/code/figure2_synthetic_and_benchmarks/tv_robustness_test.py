"""
TV-robustness test: does turning on tv_lambda destroy sequential ordering on
Figure 2 Systems 1a and 1b (ground-truth sequential handoffs)?

If matched Kendall τ stays near 1.0 across tv_lambda ∈ {0, 0.01, 0.1, 1.0, 10.0},
TV is a safe opt-in — recovers concurrent where it exists (already shown on
System 2, matched CP 0.05 → 0.49) without inducing spurious concurrency on
sequential data.

If matched τ collapses toward 0 under TV on sequential systems, TV is a bad
estimator regardless of the concurrent gain, and cannot be safely applied to
Fig 5 to test the A1→A2 handoff.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import torch
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import linear_sum_assignment
from scipy.stats import kendalltau

REPO = Path(__file__).resolve().parents[3]  # → scJDO repo root
sys.path.insert(0, str(REPO))
from scjdo.archetypes.decompose import jacobian_modes

OUT_JSON = REPO / "analysis" / "results" / "figure2" / "tv_robustness_test.json"
_RESULTS: dict = {}

METHODS_N_WINDOWS = 100
METHODS_N_PCS = 50
RANK = 5
N_SEEDS = 8
GLOBAL_SEED = 20260515
TENSOR_NOISE = 0.010


def stable_seed(*items, base=GLOBAL_SEED):
    p = "|".join(map(str, items)).encode()
    d = hashlib.sha256(p).hexdigest()[:12]
    return (base + int(d, 16)) % (2**32 - 1)


def normalize01(x, axis=0):
    x = np.asarray(x, dtype=np.float32)
    return (x - x.min(axis=axis, keepdims=True)) / (
        x.max(axis=axis, keepdims=True) - x.min(axis=axis, keepdims=True) + 1e-8
    )


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def plateau_window(t, start, end, edge=0.012):
    return sigmoid((t - start) / edge) * sigmoid((end - t) / edge)


def system1a_profiles(seed_index):
    rng = np.random.default_rng(stable_seed("system_1a_sharp_step_handoff", seed_index, "profiles"))
    t = np.linspace(0, 1, METHODS_N_WINDOWS, dtype=np.float32)
    centers = np.linspace(0.13, 0.87, RANK) + rng.normal(0, 0.006, RANK)
    widths = np.array([0.055, 0.060, 0.058, 0.060, 0.055]) * rng.uniform(0.92, 1.08, RANK)
    A = []
    for c, w in zip(centers, widths):
        prof = plateau_window(t, c - w, c + w, edge=0.008)
        prof += 0.012 * gaussian_filter1d(rng.standard_normal(METHODS_N_WINDOWS), sigma=1.0)
        A.append(prof)
    return t, normalize01(np.stack(A, axis=1), axis=0).astype(np.float32)


def system1b_profiles(seed_index):
    rng = np.random.default_rng(stable_seed("system_1b_gradual_sigmoid_handoff", seed_index, "profiles"))
    t = np.linspace(0, 1, METHODS_N_WINDOWS, dtype=np.float32)
    centers = np.linspace(0.16, 0.84, RANK) + rng.normal(0, 0.010, RANK)
    widths = np.array([0.14, 0.17, 0.19, 0.17, 0.14]) * rng.uniform(0.92, 1.10, RANK)
    A = []
    for c, w in zip(centers, widths):
        prof = plateau_window(t, c - w, c + w, edge=0.055)
        prof += 0.018 * np.sin(2 * np.pi * (t + c))
        prof += 0.010 * gaussian_filter1d(rng.standard_normal(METHODS_N_WINDOWS), sigma=2.0)
        A.append(prof)
    return t, normalize01(np.stack(A, axis=1), axis=0).astype(np.float32)


def make_operator_basis(seed, d=METHODS_N_PCS, rank=RANK):
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


def match_on_patterns(patterns_gt, patterns_rec):
    """Match by pattern (operator) cosine — independent of activation profiles."""
    K = patterns_gt.shape[0]
    g = patterns_gt.reshape(K, -1)
    r = patterns_rec.reshape(K, -1)
    g = g / (np.linalg.norm(g, axis=1, keepdims=True) + 1e-8)
    r = r / (np.linalg.norm(r, axis=1, keepdims=True) + 1e-8)
    sim = g @ r.T  # (K, K)
    row, col = linear_sum_assignment(-sim)
    return col, sim


def kendall_index_vs_peak(A, t):
    peaks = t[np.argmax(A, axis=0)].astype(float)
    tau = kendalltau(np.arange(len(peaks)), peaks).statistic
    return float(0.0 if not np.isfinite(tau) else tau)


def run(system_name, profile_fn):
    print(f"\n{'=' * 76}")
    print(f"System: {system_name}   (ground-truth sequential handoff, GT Kendall τ = 1)")
    print("=" * 76)
    print(f"{'tv_λ':>6}  {'matched τ mean':>14}  {'median':>7}  {'N(τ≥0.9)/24':>12}  {'recon_err':>10}")
    _RESULTS[system_name] = {}
    for tv in (0.0, 0.01, 0.1, 1.0, 10.0):
        taus = []
        errs = []
        for si in range(N_SEEDS):
            t, A_gt = profile_fn(si)
            J, B = make_tensor(A_gt, seed=stable_seed(system_name, si, "tensor"))
            patterns, activations, err = jacobian_modes(
                torch.from_numpy(J), rank=RANK, n_restarts=5,
                seed=stable_seed(system_name, si, "decomp_nmf"),
                tv_lambda=tv,
            )
            # Match on PATTERNS (operators B_k), not activations.
            perm, _ = match_on_patterns(B, patterns.numpy())
            A_rec = normalize01(activations.numpy(), axis=0)
            A_rec_matched = A_rec[:, perm]
            taus.append(kendall_index_vs_peak(A_rec_matched, t))
            errs.append(float(err))
        n_hit = int(sum(1 for x in taus if x >= 0.9))
        print(f"{tv:>6}  {np.mean(taus):>14.3f}  {np.median(taus):>7.3f}  "
              f"{n_hit:>7d}/{N_SEEDS:<4d}  {np.mean(errs):>10.4f}")
        _RESULTS[system_name][f"tv_lambda={tv}"] = {
            "matched_tau_mean": float(np.mean(taus)),
            "matched_tau_median": float(np.median(taus)),
            "n_seeds_ge_0.9": n_hit,
            "n_seeds": N_SEEDS,
            "recon_err_mean": float(np.mean(errs)),
            "per_seed_tau": [float(x) for x in taus],
        }


def main():
    run("system_1a_sharp_step_handoff", system1a_profiles)
    run("system_1b_gradual_sigmoid_handoff", system1b_profiles)
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(_RESULTS, f, indent=2)
    print(f"\nSaved: {OUT_JSON}")


if __name__ == "__main__":
    main()
