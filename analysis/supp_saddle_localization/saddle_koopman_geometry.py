"""
Toggle-switch pitchfork — Koopman geometry-descriptor readouts.

Extends the four-readout saddle analysis (Re(λ_max), tr(J), P⊥JP⊥, FTLE
— all reported failing to localise the analytic pitchfork at τ_crit ≈
0.004 in `saddle_readouts.py`) with four *branch-free* geometry
descriptors computed on the local Koopman operator over the same
pseudotime grid:

    R5  henrici        ν = sqrt(||K||_F² − Σ|λ|²) / ||K||_F, in [0, 1]
    R6  reactivity     σ_max(K) / ρ(K), ≥ 1, tight iff normal
    R7  transient_gain max_n ||K^n||_2 over n ≤ 30
    R8  eigvec_cond    κ(V), large when K is nearly defective

Same synthetic system as `saddle_readouts.py`, same seed, same drift
model, same τ_crit and ±0.05 tolerance. Result either:
- All four geometry descriptors also miss τ_crit → strengthens the
  manuscript's readout-agnostic saddle argument by adding an
  independent, log-free family of readouts to the failing set.
- One or more descriptors localise τ_crit where the four existing
  readouts do not → genuinely new positive result; the geometry
  axis surfaces bifurcation structure that the eigenvalue-based
  axis cannot.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path("/Users/terooatt/Downloads/scJDO")
sys.path.insert(0, str(REPO))

# Reuse the exact toggle-switch construction so the two scripts share
# the analytic ground truth and the trained drift regime.
from analysis.supp_saddle_localization.saddle_readouts import (
    ALPHA_CRIT, ALPHA_MIN, ALPHA_MAX, TAU_CRIT, TOLERANCE,
    N_CELLS, N_LATENT, N_EPOCHS, N_GRID, SEED, BANDWIDTH,
    build_adata, interior_argmax, interior_argmin,
)
from scjdo.tl import fit_drift


def score(tau_peak: float) -> dict:
    err = abs(float(tau_peak) - float(TAU_CRIT))
    return {"tau_peak": float(tau_peak), "err": float(err),
            "pass": bool(err <= TOLERANCE)}


def main() -> None:
    print(f"τ_crit (analytic) = {TAU_CRIT:.4f}   tolerance = ±{TOLERANCE}")
    print(f"scoring window   = [{max(0.0, TAU_CRIT - TOLERANCE):.3f}, "
          f"{TAU_CRIT + TOLERANCE:.3f}]")
    print()

    adata, _, _ = build_adata(SEED)
    print(f"cells = {adata.n_obs}   latent = {N_LATENT}")

    t0 = time.time()
    fit_drift(
        adata,
        rep="X_pca",
        time_key="pseudotime",
        n_epochs=N_EPOCHS,
        n_archetypes=4,
        n_eff_min=20.0,
        n_boot=10,
        grid_size=N_GRID,
        seed=SEED,
        archetype_method="koopman",
        koopman_kwargs=dict(mode="local", window_half=8, ridge=1e-4,
                            whiten=False),
        key_added="scjdo_koop",
        verbose=False,
    )
    print(f"fit_drift(koopman) in {time.time() - t0:.1f}s")

    res = adata.uns["scjdo_koop"]
    grid = np.asarray(res["t_centers"], dtype=np.float64)
    koop = res["koopman"]
    g = koop["geometry"]
    dt = np.asarray(koop["delta_tau"])

    # Sanity — same drift field must reproduce the R1 baseline failure
    # (peak deep in committed-branch region, not at τ_crit ≈ 0.004).
    lam_max = np.asarray(res["max_real_eig"], dtype=np.float64)
    tau_R1 = interior_argmax(lam_max, grid)
    print(f"R1  Re(λ_max) baseline peak τ = {tau_R1:.3f}    "
          f"err = {abs(tau_R1 - TAU_CRIT):.3f}  "
          f"({'PASS' if abs(tau_R1 - TAU_CRIT) <= TOLERANCE else 'FAIL'})")
    print()

    # Score each geometry descriptor — argmax on the interior window,
    # same convention the four existing readouts use.
    def _peak(field: str) -> float:
        arr = np.asarray(g[field], dtype=np.float64)
        return interior_argmax(arr, grid)

    tau_R5 = _peak("henrici")
    tau_R6 = _peak("reactivity")
    tau_R7 = _peak("transient_gain")
    tau_R8 = _peak("eigvec_cond")

    results = {
        "tau_crit":       float(TAU_CRIT),
        "tolerance":      float(TOLERANCE),
        "R1_lambda_max":  score(tau_R1),
        "R5_henrici":     score(tau_R5),
        "R6_reactivity":  score(tau_R6),
        "R7_transient_gain": score(tau_R7),
        "R8_eigvec_cond": score(tau_R8),
        "delta_tau_median": float(np.median(dt[dt > 0])),
        "nyquist_ang_median": float(np.median(np.pi / dt[dt > 0])),
        "notes": (
            "Geometry descriptors are branch-free (no log λ), so they "
            "sidestep the Nyquist ambiguity that constrains R1/R3-style "
            "eigenvalue readouts. Result is directly comparable to "
            "saddle_localization_v3_metrics.json / saddle_readouts.json."
        ),
    }

    print("=" * 72)
    print(f"{'Readout':32s} {'τ_peak':>8s} {'err vs τ_crit':>14s} {'PASS?':>6s}")
    print("-" * 72)
    for k in ["R1_lambda_max", "R5_henrici", "R6_reactivity",
              "R7_transient_gain", "R8_eigvec_cond"]:
        r = results[k]
        print(f"{k:32s} {r['tau_peak']:>8.3f} {r['err']:>14.3f} "
              f"{'PASS' if r['pass'] else 'FAIL':>6s}")
    print("=" * 72)

    # ── Verdict summary ──────────────────────────────────────────────────
    geom_keys = ["R5_henrici", "R6_reactivity", "R7_transient_gain", "R8_eigvec_cond"]
    n_pass = sum(1 for k in geom_keys if results[k]["pass"])
    print()
    if n_pass == 0:
        print("VERDICT — All four Koopman geometry descriptors also miss τ_crit.")
        print("         → Extends the readout-agnostic saddle limit; the failure")
        print("           is upstream of every readout tested (R1..R8), not")
        print("           specific to eigenvalue-based summaries.")
    elif n_pass == len(geom_keys):
        print("VERDICT — All four Koopman geometry descriptors localise τ_crit.")
        print("         → Positive result: the geometry axis recovers saddle")
        print("           structure the eigenvalue axis cannot.")
    else:
        which = [k for k in geom_keys if results[k]["pass"]]
        print(f"VERDICT — {n_pass}/{len(geom_keys)} Koopman geometry descriptors "
              f"localise τ_crit: {which}")
        print("         → Partial positive: report the specific descriptor(s)")
        print("           that recover the pitchfork as a candidate saddle readout.")

    outdir = Path("/Users/terooatt/Downloads/scJDO/scratchpad_run")
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "saddle_koopman_geometry.json").write_text(json.dumps(results, indent=2))

    # Save the curves for optional plotting.
    np.savez_compressed(
        outdir / "saddle_koopman_geometry_curves.npz",
        grid=grid,
        henrici=np.asarray(g["henrici"]),
        reactivity=np.asarray(g["reactivity"]),
        transient_gain=np.asarray(g["transient_gain"]),
        eigvec_cond=np.asarray(g["eigvec_cond"]),
        R1_lambda_max=lam_max,
        delta_tau=dt,
    )
    print(f"\nSaved: {outdir / 'saddle_koopman_geometry.json'}")
    print(f"Saved: {outdir / 'saddle_koopman_geometry_curves.npz'}")


if __name__ == "__main__":
    main()
