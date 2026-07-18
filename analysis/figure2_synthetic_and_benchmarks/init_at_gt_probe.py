"""
Test #5 from critique: is the concurrent-recovery failure an optimizer
failure or an objective preference?

Method: initialize the semi-NMF W at the ground-truth A_gt on Fig-2-System-2
seeds and let it iterate. Compare final residual to (a) the residual at GT
itself (0.1103), and (b) the residual of a cold-started run (0.1105).

  If W-init-at-GT converges to GT residual (0.1103): objective is near-flat,
    and cold starts drift to a wider basin around spread-peaks. → basin
    issue, more restarts / smarter init would fix; "objective indifferent"
    wording is fine as a first-order description.

  If W-init-at-GT drifts away from GT to something worse: the objective
    actively pushes W away from the concurrent solution — that is not
    indifference. → decomposition-choice-property language stands but the
    residual-ratio argument does not.

  If W-init-at-GT drifts to 0.1105 exactly: the cold-start and warm-start
    are hitting the same optimum. Then "spread solution has lower residual
    than GT" is the real finding, and my earlier text is backwards.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import torch
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import nnls

REPO = Path("/Users/terooatt/Downloads/scJDO")
sys.path.insert(0, str(REPO))
from scjdo.archetypes.decompose import _nnls_rows, _rescale

OUT_JSON = REPO / "analysis" / "results" / "figure2" / "init_at_gt_probe.json"

METHODS_N_WINDOWS = 100
METHODS_N_PCS = 50
RANK_TRUE = 5
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


def gaussian_bump(t, center, width):
    return np.exp(-0.5 * ((t - center) / width) ** 2)


def system2_profiles(seed_index):
    rng = np.random.default_rng(stable_seed("system_2_concurrent_varying_overlap", seed_index, "profiles"))
    t = np.linspace(0, 1, METHODS_N_WINDOWS, dtype=np.float32)
    center = 0.52 + rng.normal(0, 0.006)
    widths = np.array([0.075, 0.105, 0.140, 0.180, 0.220]) * rng.uniform(0.88, 1.12, RANK_TRUE)
    A = []
    for w in widths:
        local_center = center + rng.normal(0, 0.008)
        prof = gaussian_bump(t, local_center, w)
        prof += 0.020 * gaussian_bump(t, local_center + rng.normal(0, 0.006), w * 0.60)
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


def concurrent_prevalence(A, t, gap=0.045):
    peaks = t[np.argmax(A, axis=0)]
    K = A.shape[1]
    if K < 2:
        return 0.0
    hit = pairs = 0
    for i in range(K):
        for j in range(i+1, K):
            pairs += 1
            if abs(peaks[i] - peaks[j]) <= gap:
                hit += 1
    return hit / pairs


def als_from_init(M, W_init, max_iter=500, tol=1e-6):
    """Run the exact same ALS as _seminmf but with a user-provided W_init."""
    T, D = M.shape
    W = W_init.copy().astype(np.float32)
    rank = W.shape[1]
    prev_err = np.inf
    err_trace = []
    for it in range(max_iter):
        WtW = W.T @ W + 1e-6 * np.eye(rank)
        H = np.linalg.solve(WtW, W.T @ M)
        W = _nnls_rows(H, M)
        err = float(np.linalg.norm(M - W @ H, "fro"))
        err_trace.append(err)
        if abs(prev_err - err) / (prev_err + 1e-8) < tol:
            break
        prev_err = err
    W, H = _rescale(W, H)
    return W, H, err, err_trace


def main():
    print("Init-at-GT probe on Fig-2-System-2")
    print("=" * 88)
    print(f"{'seed':>4}  {'residual @ W=A_gt (LS-H)':>26}  {'residual @ cold-start':>22}  "
          f"{'residual @ init=A_gt (ALS)':>28}  {'cp @ init=A_gt':>16}")
    print("-" * 88)
    r_gt_ls, r_cold, r_warm, cp_warm = [], [], [], []
    for si in range(N_SEEDS):
        t, A_gt = system2_profiles(si)
        J, B = make_tensor(A_gt, seed=stable_seed("system_2_concurrent_varying_overlap", si, "tensor"))
        T = A_gt.shape[0]
        d = J.shape[1]
        M = J.reshape(T, d*d).astype(np.float32)

        # (1) Static: residual with W fixed at A_gt and H at LS-best given W.
        H_from_gt = np.linalg.pinv(A_gt) @ M
        res_gt_ls = float(np.linalg.norm(M - A_gt @ H_from_gt, "fro"))

        # (2) Cold-start ALS (single restart for comparability).
        rng = np.random.default_rng(stable_seed("system_2_concurrent_varying_overlap", si, "cold"))
        W_cold = np.abs(rng.standard_normal((T, RANK_TRUE)).astype(np.float32)) + 1e-4
        W_cold /= W_cold.sum(axis=0, keepdims=True) + 1e-8
        _, _, res_cold, _ = als_from_init(M, W_cold)

        # (3) Warm-start ALS: initialise W at the ground-truth A_gt.
        W_warm, H_warm, res_warm, trace = als_from_init(M, A_gt.copy())
        A_after = normalize01(W_warm, axis=0)
        cp = concurrent_prevalence(A_after, t)

        r_gt_ls.append(res_gt_ls)
        r_cold.append(res_cold)
        r_warm.append(res_warm)
        cp_warm.append(cp)

        print(f"{si:>4}  {res_gt_ls:>26.4f}  {res_cold:>22.4f}  "
              f"{res_warm:>28.4f}  {cp:>16.3f}")

    print("=" * 88)
    print(f"MEAN  {np.mean(r_gt_ls):>26.4f}  {np.mean(r_cold):>22.4f}  "
          f"{np.mean(r_warm):>28.4f}  {np.mean(cp_warm):>16.3f}")
    print()
    print("Interpretation:")
    if np.mean(r_warm) > np.mean(r_gt_ls) + 1e-4:
        print("  Warm-started ALS drifts AWAY from GT residual → the objective is not merely")
        print("  indifferent; it actively prefers a solution worse than GT. This is decomposition-")
        print("  choice behaviour (algorithm-level) but not 'flat objective'.")
    elif np.mean(r_warm) <= np.mean(r_gt_ls) + 1e-4:
        print("  Warm-started ALS stays at (or improves on) GT residual → cold starts land in a")
        print("  wider basin around spread-peaks. 'Near-flat with large spread-basin' framing is")
        print("  supported. Adding a temporal prior tips the objective toward GT.")

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "n_seeds": N_SEEDS,
        "system": "system_2_concurrent_varying_overlap",
        "residual_at_W_eq_A_gt_LS": {
            "per_seed": [float(x) for x in r_gt_ls],
            "mean": float(np.mean(r_gt_ls)),
        },
        "residual_cold_start_ALS": {
            "per_seed": [float(x) for x in r_cold],
            "mean": float(np.mean(r_cold)),
        },
        "residual_warm_start_ALS_init_at_A_gt": {
            "per_seed": [float(x) for x in r_warm],
            "mean": float(np.mean(r_warm)),
        },
        "concurrent_prevalence_warm_start": {
            "per_seed": [float(x) for x in cp_warm],
            "mean": float(np.mean(cp_warm)),
        },
    }
    with open(OUT_JSON, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"\nSaved: {OUT_JSON}")


if __name__ == "__main__":
    main()
