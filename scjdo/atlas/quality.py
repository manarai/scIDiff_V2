"""
Drift-field quality gate for scOpAtlas.

Task 3.3. Before this module, StableOperatorAtlas classified ANY drift
field however badly fit. The classifier always returns a partition, so a
field that carries an arbitrary rotational component still produces
clean-looking regime maps on a UMAP — the failure is silent, and silent
failures that produce publishable-looking figures are the worst category.

The gate lives here (a diagnostic property of the trained drift field)
but is enforced at ``StableOperatorAtlas.build`` — a gate that lives in
one place and is checked in another is a gate people route around.

Three criteria, each SKIPPED gracefully when its inputs aren't
available, PASSED/FAILED individually when they are:

1. **Transport error (pseudotime holdout, no refit).** Bin cells by
   pseudotime into K bins. For each consecutive pair, push bin_i cells
   forward by Δt under the drift and compare the pushed distribution to
   the actual bin_{i+1} distribution via sliced Wasserstein distance.
   PASS iff the model's pushed distribution is closer than the no-op
   "stay put" baseline. The work order asks for refit-on-subset, which
   isn't feasible inside a check that receives an already-trained model;
   this pseudotime-holdout variant captures the same spirit without
   coupling the check to a specific trainer.

2. **Velocity-prior agreement.** When adata carries an RNA-velocity
   estimate in the same latent space, mean per-cell cosine similarity
   between drift(x, t) and velocity. Never expected to be high — RNA
   velocity is noisy and often only weakly aligned. Threshold is loose;
   the point is "non-random," not "well-aligned."

3. **Ensemble regime concordance.** Uses the OperatorEnsemble output
   (`adata.obs['regime_concordance']`) — the most direct measure of
   gauge robustness. PASS iff the discordant fraction (cells below the
   concordance threshold) is below the tolerated fraction.

The overall gate PASSES iff at least one criterion actually ran AND every
criterion that ran passed. On failure, ``StableOperatorAtlas.build`` raises
by default; pass ``force=True`` to override with a loud banner. The full
report is always written to ``adata.uns['drift_quality']`` so it travels
with the object.
"""

from __future__ import annotations

from typing import Dict, Optional

import numpy as np
import torch


# ---------------------------------------------------------------------------
# Sliced Wasserstein distance (pure numpy — no ot library dependency).
# ---------------------------------------------------------------------------

def sliced_wasserstein(
    X1: np.ndarray,
    X2: np.ndarray,
    n_projections: int = 50,
    rng: Optional[np.random.Generator] = None,
) -> float:
    """1D-averaged Wasserstein distance across random unit projections.

    O(n_projections * (n1 + n2) log(n1 + n2)). Fine for the sizes this
    module deals with (per-bin ~ hundreds of cells).
    """
    if rng is None:
        rng = np.random.default_rng(0)
    X1 = np.asarray(X1)
    X2 = np.asarray(X2)
    if X1.size == 0 or X2.size == 0:
        return float("nan")
    d = X1.shape[1]
    proj = rng.standard_normal((n_projections, d))
    proj /= np.linalg.norm(proj, axis=1, keepdims=True) + 1e-12

    total = 0.0
    for p in proj:
        u = np.sort(X1 @ p)
        v = np.sort(X2 @ p)
        # Quantile match via linear interpolation onto a common grid.
        n = min(len(u), len(v))
        qs = np.linspace(0, 1, n)
        uq = np.quantile(u, qs)
        vq = np.quantile(v, qs)
        total += float(np.mean(np.abs(uq - vq)))
    return total / n_projections


# ---------------------------------------------------------------------------
# Criterion 1: transport error (pseudotime-holdout, no refit)
# ---------------------------------------------------------------------------

def _push_forward(model, x: torch.Tensor, t_start: float, dt: float) -> np.ndarray:
    """One Euler step: x + f(x, t_start) * dt.

    A single step is enough for the diagnostic — the point is to check
    the drift's DIRECTION and magnitude are sensible in latent space,
    not to build a full integrator inside the quality gate.
    """
    with torch.no_grad():
        t = torch.full((x.shape[0],), float(t_start), dtype=x.dtype, device=x.device)
        f = model(x, t)
    return (x + dt * f).cpu().numpy()


def transport_holdout_score(
    model,
    X: np.ndarray,
    pseudotime: np.ndarray,
    n_bins: int = 5,
    n_projections: int = 50,
    device: str = "cpu",
    rng: Optional[np.random.Generator] = None,
) -> Dict[str, float]:
    """Predict-vs-observed sliced-Wasserstein by pseudotime bin.

    Returns dict:
        model_swd  — mean SWD(pushed bin_i, actual bin_{i+1}) across pairs.
        baseline_swd — mean SWD(bin_i, actual bin_{i+1}) across pairs
            (no drift applied). If the model's pushed distribution is
            farther than this, the drift is doing harm relative to
            "stay put".
        ratio — model_swd / baseline_swd. < 1.0 means the drift is
            informative for the transport direction.
        n_pairs — number of consecutive-bin pairs used.
    """
    if rng is None:
        rng = np.random.default_rng(0)

    order = np.argsort(pseudotime)
    Xo = np.asarray(X)[order]
    to = np.asarray(pseudotime)[order]
    n = len(Xo)
    bin_size = max(2, n // n_bins)

    bins = []
    for i in range(n_bins):
        lo, hi = i * bin_size, (i + 1) * bin_size if i < n_bins - 1 else n
        if hi - lo >= 2:
            bins.append((Xo[lo:hi], to[lo:hi]))

    if len(bins) < 2:
        return dict(model_swd=float("nan"), baseline_swd=float("nan"),
                    ratio=float("nan"), n_pairs=0)

    model_swds, baseline_swds = [], []
    for i in range(len(bins) - 1):
        X_i, t_i = bins[i]
        X_j, _ = bins[i + 1]
        t_mid = float(np.mean(t_i))
        dt = float(np.mean(bins[i + 1][1]) - t_mid)
        if dt <= 0:
            continue

        x_tensor = torch.tensor(X_i, dtype=torch.float32, device=device)
        pushed = _push_forward(model, x_tensor, t_mid, dt)

        model_swds.append(sliced_wasserstein(pushed, X_j, n_projections, rng))
        baseline_swds.append(sliced_wasserstein(X_i, X_j, n_projections, rng))

    if not model_swds:
        return dict(model_swd=float("nan"), baseline_swd=float("nan"),
                    ratio=float("nan"), n_pairs=0)

    m = float(np.mean(model_swds))
    b = float(np.mean(baseline_swds))
    return dict(
        model_swd=m,
        baseline_swd=b,
        ratio=(m / b if b > 0 else float("inf")),
        n_pairs=len(model_swds),
    )


# ---------------------------------------------------------------------------
# Criterion 2: velocity-prior agreement
# ---------------------------------------------------------------------------

def _find_velocity(adata) -> Optional[np.ndarray]:
    """Locate an RNA-velocity array in the same latent space as X.

    Fallback chain:
      - adata.obsm['velocity']         (same coords as X)
      - adata.obsm['velocity_pca']     (scvelo convention on PC space)
      - None (skip the check)
    """
    for key in ("velocity", "velocity_pca"):
        if key in adata.obsm:
            v = np.asarray(adata.obsm[key])
            return v
    return None


def velocity_agreement_score(
    model,
    X: np.ndarray,
    pseudotime: np.ndarray,
    velocity: np.ndarray,
    device: str = "cpu",
    subsample: Optional[int] = 2000,
    rng: Optional[np.random.Generator] = None,
) -> Dict[str, float]:
    """Mean per-cell cosine similarity between drift and velocity.

    Returns dict:
        mean_cos — mean cosine over cells with defined velocity.
        median_cos, frac_positive, n_used.
    """
    if rng is None:
        rng = np.random.default_rng(0)

    X = np.asarray(X)
    velocity = np.asarray(velocity)
    if X.shape != velocity.shape:
        raise ValueError(
            f"velocity shape {velocity.shape} != X shape {X.shape}; "
            "project velocity to the same latent space as X first."
        )

    # Drop cells with all-zero or all-NaN velocity.
    finite = np.isfinite(velocity).all(axis=1)
    nonzero = np.linalg.norm(velocity, axis=1) > 1e-12
    valid = finite & nonzero
    if valid.sum() == 0:
        return dict(mean_cos=float("nan"), median_cos=float("nan"),
                    frac_positive=float("nan"), n_used=0)

    idx = np.flatnonzero(valid)
    if subsample is not None and len(idx) > subsample:
        idx = rng.choice(idx, size=subsample, replace=False)

    X_sub = torch.tensor(X[idx], dtype=torch.float32, device=device)
    t_sub = torch.tensor(pseudotime[idx], dtype=torch.float32, device=device)
    with torch.no_grad():
        drift = model(X_sub, t_sub).cpu().numpy()

    v_sub = velocity[idx]
    dot = np.sum(drift * v_sub, axis=1)
    norms = np.linalg.norm(drift, axis=1) * np.linalg.norm(v_sub, axis=1)
    cos = np.where(norms > 1e-12, dot / np.maximum(norms, 1e-12), 0.0)
    return dict(
        mean_cos=float(cos.mean()),
        median_cos=float(np.median(cos)),
        frac_positive=float((cos > 0).mean()),
        n_used=int(len(cos)),
    )


# ---------------------------------------------------------------------------
# Criterion 3: ensemble concordance
# ---------------------------------------------------------------------------

def ensemble_concordance_score(
    concordance: np.ndarray,
    threshold: float = 0.5,
) -> Dict[str, float]:
    """Discordant fraction from precomputed concordance."""
    concordance = np.asarray(concordance, dtype=float)
    if concordance.size == 0:
        return dict(discordant_fraction=float("nan"),
                    mean_concordance=float("nan"), n_cells=0)
    return dict(
        discordant_fraction=float((concordance < threshold).mean()),
        mean_concordance=float(concordance.mean()),
        n_cells=int(len(concordance)),
    )


# ---------------------------------------------------------------------------
# Top-level entry
# ---------------------------------------------------------------------------

def check_drift_quality(
    model,
    adata,
    use_rep: str = "X_pca",
    pseudotime_key: str = "pseudotime",
    concordance_key: str = "regime_concordance",
    n_bins: int = 5,
    n_projections: int = 50,
    subsample_velocity: Optional[int] = 2000,
    transport_ratio_max: float = 0.95,
    velocity_cos_min: float = 0.05,
    concordance_threshold: float = 0.5,
    max_discordant_fraction: float = 1.0 / 3.0,
    device: str = "cpu",
    seed: int = 0,
) -> Dict:
    """
    Run the three-criterion drift-quality gate. Returns a report dict
    with an overall ``passed`` field plus per-criterion detail.

    Pass conditions:
      - Transport error: model_swd / baseline_swd < transport_ratio_max
        (default 0.95 — model beats no-op stay-put by 5%). Skipped if
        insufficient bins.
      - Velocity agreement: mean cosine > velocity_cos_min (default 0.05
        — weak-but-positive alignment). Skipped if no velocity in adata.
      - Ensemble concordance: discordant fraction < max_discordant_fraction
        (default 1/3). Skipped if regime_concordance is not in
        adata.obs (i.e. no ensemble was run).

    Overall PASS iff at least one criterion ran AND every criterion that
    ran passed.
    """
    rng = np.random.default_rng(seed)
    report: Dict = {"criteria": {}, "thresholds": dict(
        transport_ratio_max=transport_ratio_max,
        velocity_cos_min=velocity_cos_min,
        concordance_threshold=concordance_threshold,
        max_discordant_fraction=max_discordant_fraction,
    )}

    if use_rep not in adata.obsm:
        raise ValueError(f"'{use_rep}' not in adata.obsm")
    if pseudotime_key not in adata.obs.columns:
        raise ValueError(f"'{pseudotime_key}' not in adata.obs")
    X = np.asarray(adata.obsm[use_rep])
    pseudotime = np.asarray(adata.obs[pseudotime_key].values, dtype=float)

    # --- 1. transport ---
    t_score = transport_holdout_score(
        model, X, pseudotime, n_bins=n_bins,
        n_projections=n_projections, device=device, rng=rng,
    )
    if t_score["n_pairs"] > 0:
        ran = True
        passed = t_score["ratio"] < transport_ratio_max
        note = (f"model SWD {t_score['model_swd']:.4f} vs baseline "
                f"{t_score['baseline_swd']:.4f} (ratio {t_score['ratio']:.3f}; "
                f"pass < {transport_ratio_max})")
    else:
        ran = False; passed = False
        note = "not enough pseudotime bins with cells; check skipped."
    report["criteria"]["transport_error"] = dict(
        ran=ran, passed=passed, note=note, **t_score,
    )

    # --- 2. velocity ---
    velocity = _find_velocity(adata)
    if velocity is None:
        report["criteria"]["velocity_agreement"] = dict(
            ran=False, passed=False,
            note="no adata.obsm['velocity'] or 'velocity_pca'; check skipped.",
        )
    else:
        v_score = velocity_agreement_score(
            model, X, pseudotime, velocity,
            device=device, subsample=subsample_velocity, rng=rng,
        )
        passed = np.isfinite(v_score["mean_cos"]) and v_score["mean_cos"] > velocity_cos_min
        note = (f"mean cos = {v_score['mean_cos']:.3f} over {v_score['n_used']} "
                f"cells (pass > {velocity_cos_min})")
        report["criteria"]["velocity_agreement"] = dict(
            ran=True, passed=bool(passed), note=note, **v_score,
        )

    # --- 3. ensemble concordance ---
    if concordance_key not in adata.obs.columns:
        report["criteria"]["ensemble_concordance"] = dict(
            ran=False, passed=False,
            note=(f"no adata.obs['{concordance_key}']; run OperatorEnsemble.build "
                  "first to enable this check."),
        )
    else:
        e_score = ensemble_concordance_score(
            adata.obs[concordance_key].values,
            threshold=concordance_threshold,
        )
        passed = np.isfinite(e_score["discordant_fraction"]) and \
                 e_score["discordant_fraction"] < max_discordant_fraction
        note = (f"discordant fraction = {e_score['discordant_fraction']:.2%} "
                f"(pass < {max_discordant_fraction:.2%})")
        report["criteria"]["ensemble_concordance"] = dict(
            ran=True, passed=bool(passed), note=note, **e_score,
        )

    # --- overall ---
    ran_criteria = [c for c in report["criteria"].values() if c["ran"]]
    n_ran = len(ran_criteria)
    passed_all = all(c["passed"] for c in ran_criteria) if ran_criteria else False
    report["n_criteria_ran"] = n_ran
    report["passed"] = bool(n_ran > 0 and passed_all)
    report["summary"] = _summary_string(report)
    return report


def _summary_string(report: Dict) -> str:
    lines = [f"drift_quality: passed={report['passed']}  "
             f"(criteria that ran: {report['n_criteria_ran']}/3)"]
    for name, c in report["criteria"].items():
        marker = "✓" if c["passed"] else ("✗" if c["ran"] else "-")
        lines.append(f"  {marker} {name}: {c['note']}")
    return "\n".join(lines)


class DriftQualityError(RuntimeError):
    """Raised when the drift-quality gate fails and force is not set."""


def raise_if_failed(report: Dict, force: bool = False) -> None:
    """Enforcement hook — used by StableOperatorAtlas.build."""
    if report["passed"]:
        return
    if force:
        banner = (
            "\n" + "!" * 78 + "\n"
            "!! drift_quality gate FAILED but force=True — proceeding anyway.\n"
            "!! Report:\n" + report["summary"] + "\n"
            + "!" * 78 + "\n"
        )
        print(banner)
        return
    raise DriftQualityError(
        "drift_quality gate failed. Report:\n" + report["summary"]
        + "\n\nPass force=True to StableOperatorAtlas.build(...) to override "
          "(a loud banner will be printed). The full report is available "
          "in adata.uns['drift_quality'] regardless."
    )
