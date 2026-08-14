"""
scjdo.validation.identifiability
=====================================
Utilities for quantifying what is and is not identifiable in the
scJDO operator framework.

Background
----------
scJDO infers time-varying Jacobians as derivatives of a learned
drift field.  Because the drift field is estimated from static
snapshot data, the Jacobians are *not* uniquely identifiable in an
absolute sense: many different drift fields can produce similar
observed distributions but different Jacobians.

However, certain *relative* and *structural* properties of the
inferred operator space are empirically robust across model choices,
embeddings, and seeds.  This module provides tools to quantify those
properties and to test how sensitive conclusions are to model-class
variation (architecture, regularization, noise level).

Invariant properties (empirically, under prior variation within a
declared family):
    - Relative reuse of operator regimes (archetype activation ordering)
    - Low-rank structure of Jacobian evolution (variance explained by K)
    - Gene-direction / regime structure of archetypes

Non-invariant properties (explicitly acknowledged):
    - Absolute eigenvalues
    - Exact Jacobian matrix entries
    - Gene-level causal interpretation
    - Conclusions under nonlinear reparameterizations of latent space
    - Absolute timing of instability peaks (magnitude and pseudotime
      location depend on σ, spectral-norm, and velocity-guidance weight;
      previously listed as invariant, retracted per v3 benchmark showing
      mean peak-timing error 0.58 vs pre-committed ≤0.10)

Usage example
-------------
>>> from scjdo.validation.identifiability import (
...     archetype_cosine_similarity,
...     instability_peak_overlap,
...     model_sensitivity_report,
... )
>>> sim = archetype_cosine_similarity(archetypes_run1, archetypes_run2)
>>> print(f"Median cosine similarity: {sim['median']:.3f}")
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import torch

# ---------------------------------------------------------------------------
# Archetype identity across runs / model choices
# ---------------------------------------------------------------------------

def archetype_cosine_similarity(
    archetypes_a: torch.Tensor,
    archetypes_b: torch.Tensor,
) -> dict:
    """
    Compute cosine similarity between matched archetype pairs from two runs.

    Archetypes are matched greedily by maximum cosine similarity (Hungarian
    matching is not required because K is small and similarity is typically
    unambiguous).

    Parameters
    ----------
    archetypes_a : Tensor, shape (K, D, D)
        Archetype matrices from run A.
    archetypes_b : Tensor, shape (K, D, D)
        Archetype matrices from run B.

    Returns
    -------
    dict with keys:
        per_archetype : list[float]  cosine similarity for each matched pair
        median        : float        median across all pairs
        min           : float        minimum (worst-case pair)
        summary       : str
    """
    Ka = archetypes_a.shape[0]
    Kb = archetypes_b.shape[0]
    K = min(Ka, Kb)

    # Flatten to vectors
    flat_a = archetypes_a[:Ka].reshape(Ka, -1).float()
    flat_b = archetypes_b[:Kb].reshape(Kb, -1).float()

    # Normalise
    norm_a = flat_a / (flat_a.norm(dim=1, keepdim=True) + 1e-8)
    norm_b = flat_b / (flat_b.norm(dim=1, keepdim=True) + 1e-8)

    # Similarity matrix (Ka x Kb)
    sim_matrix = (norm_a @ norm_b.T).numpy()

    # Greedy matching
    matched: list[float] = []
    used_b: set[int] = set()
    for i in range(K):
        row = sim_matrix[i].copy()
        for j in used_b:
            row[j] = -np.inf
        best_j = int(np.argmax(row))
        matched.append(float(sim_matrix[i, best_j]))
        used_b.add(best_j)

    median_sim = float(np.median(matched))
    min_sim = float(np.min(matched))

    summary = (
        f"Archetype cosine similarity (K={K})\n"
        f"  Per-archetype: {[f'{v:.3f}' for v in matched]}\n"
        f"  Median: {median_sim:.3f}  Min: {min_sim:.3f}"
    )

    return {
        "per_archetype": matched,
        "median": median_sim,
        "min": min_sim,
        "summary": summary,
    }


# ---------------------------------------------------------------------------
# Instability peak timing across runs
# ---------------------------------------------------------------------------

def instability_peak_overlap(
    instability_curves: list[np.ndarray],
    peak_window: float = 0.1,
    ground_truth_peak: Optional[float] = None,
) -> dict:
    """
    Measure how consistently instability peaks are located across runs.

    Cross-run agreement (overlap_frac) measures precision, not correctness:
    runs that agree on the wrong peak score 1.0. Pass ``ground_truth_peak``
    (e.g., the true saddle crossing on a synthetic benchmark) to also
    report per-run absolute error and the fraction of runs within
    ``peak_window`` of ground truth.

    Parameters
    ----------
    instability_curves : list of 1-D arrays
        Each array is a normalised instability curve over pseudotime [0, 1].
        All arrays should have the same length T.
    peak_window : float
        Two peaks (or a peak and ground truth) are considered "overlapping"
        if they are within this fraction of pseudotime of each other.
    ground_truth_peak : float, optional
        Known pseudotime location of the target instability event. When
        provided, the report includes accuracy against this reference in
        addition to cross-run precision.

    Returns
    -------
    dict with keys:
        peak_locations : list[float]  pseudotime of peak per run
        peak_std       : float        std of peak locations
        overlap_frac   : float        fraction of run pairs whose peaks overlap
        gt_errors      : list[float]  |peak - ground_truth| per run  (if GT given)
        gt_mean_error  : float        mean |peak - ground_truth|     (if GT given)
        gt_accuracy    : float        frac. of runs within peak_window of GT (if GT)
        summary        : str
    """
    peak_locs = []
    for curve in instability_curves:
        t = np.linspace(0, 1, len(curve))
        peak_locs.append(float(t[np.argmax(curve)]))

    peak_std = float(np.std(peak_locs))

    n = len(peak_locs)
    n_pairs = n * (n - 1) // 2
    overlap_count = 0
    for i in range(n):
        for j in range(i + 1, n):
            if abs(peak_locs[i] - peak_locs[j]) <= peak_window:
                overlap_count += 1
    overlap_frac = overlap_count / max(n_pairs, 1)

    out: dict = {
        "peak_locations": peak_locs,
        "peak_std": peak_std,
        "overlap_frac": overlap_frac,
    }

    summary_lines = [
        f"Instability peak timing ({len(instability_curves)} runs)",
        f"  Peak locations: {[f'{p:.3f}' for p in peak_locs]}",
        f"  Std of peak locations: {peak_std:.4f}",
        f"  Cross-run precision (pairs within {peak_window}): {overlap_frac:.3f}",
    ]

    if ground_truth_peak is not None:
        gt_errors = [abs(p - ground_truth_peak) for p in peak_locs]
        gt_mean_error = float(np.mean(gt_errors))
        gt_accuracy = float(np.mean([e <= peak_window for e in gt_errors]))
        out.update(
            {
                "gt_errors": gt_errors,
                "gt_mean_error": gt_mean_error,
                "gt_accuracy": gt_accuracy,
            }
        )
        summary_lines.extend(
            [
                f"  Ground-truth peak: {ground_truth_peak:.3f}",
                f"  Per-run |Δt|: {[f'{e:.3f}' for e in gt_errors]}",
                f"  Mean |Δt|: {gt_mean_error:.4f}  "
                f"Frac. within {peak_window}: {gt_accuracy:.3f}",
            ]
        )

    out["summary"] = "\n".join(summary_lines)
    return out


# ---------------------------------------------------------------------------
# Model-sensitivity report
# ---------------------------------------------------------------------------

def model_sensitivity_report(
    results_by_config: dict[str, dict],
) -> str:
    """
    Produce a formatted sensitivity report across model configurations.

    This measures ESTIMATOR PRECISION across model-class choices (depth,
    width). It does NOT measure identifiability: two architectures agreeing
    tells you the estimator is stable, not that the density constrained the
    drift. For identifiability, use ``prior_sensitivity_report`` which
    varies the actual prior knobs (σ, spectral-norm, velocity guidance,
    smoothness).

    Parameters
    ----------
    results_by_config : dict[str, dict]
        Mapping from configuration label (e.g., "depth=4, hidden=256") to
        a dict containing at minimum:
            - "archetypes"         : Tensor (K, D, D)
            - "instability_curve"  : np.ndarray (T,)
            - "auroc"              : float (optional)

    Returns
    -------
    str
        Formatted report comparing archetype similarity and instability
        peak timing across all configurations.
    """
    configs = list(results_by_config.keys())
    lines = [
        "Model sensitivity report",
        "=" * 60,
        "Invariant properties tested:",
        "  1. Archetype cosine similarity across configurations",
        "  2. Instability peak timing across configurations",
        "",
        "Non-invariant (explicitly not claimed):",
        "  - Absolute eigenvalues",
        "  - Exact Jacobian entries",
        "  - Gene-level causal interpretation",
        "=" * 60,
    ]

    # Archetype similarity: compare each config to the first
    ref_label = configs[0]
    ref_arch = results_by_config[ref_label]["archetypes"]
    lines.append(f"\nArchetype similarity vs reference ({ref_label}):")
    for label in configs[1:]:
        arch = results_by_config[label]["archetypes"]
        sim = archetype_cosine_similarity(ref_arch, arch)
        lines.append(f"  {label}: median cosine = {sim['median']:.3f}, min = {sim['min']:.3f}")

    # Instability peak timing
    curves = [results_by_config[c]["instability_curve"] for c in configs]
    peak_info = instability_peak_overlap(curves)
    lines.append("\nInstability peak timing:")
    for label, loc in zip(configs, peak_info["peak_locations"]):
        lines.append(f"  {label}: peak at pseudotime = {loc:.3f}")
    lines.append(f"  Std across configs: {peak_info['peak_std']:.4f}")
    lines.append(f"  Fraction of pairs within 0.10: {peak_info['overlap_frac']:.3f}")

    # AUROC if available
    aurocs = {
        c: results_by_config[c]["auroc"]
        for c in configs
        if "auroc" in results_by_config[c]
    }
    if aurocs:
        lines.append("\nAUROC across configurations:")
        for label, auroc in aurocs.items():
            lines.append(f"  {label}: {auroc:.4f}")
        vals = list(aurocs.values())
        lines.append(
            f"  Range: {min(vals):.4f} – {max(vals):.4f}  "
            f"(ΔAUROC = {max(vals) - min(vals):.4f})"
        )

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Prior-sensitivity: the real identifiability probe
# ---------------------------------------------------------------------------

# Canonical prior-knob families for the drift-field inverse problem.
# Sweeping architecture (depth, width, seed) at fixed prior gives you
# ESTIMATOR PRECISION. Sweeping *these* — the terms that shape which
# drift field the density selects among — gives you IDENTIFIABILITY.
PRIOR_KNOBS: tuple[str, ...] = (
    "sigma",                     # bridge/score noise level
    "spectral_norm",             # bool: Lipschitz cap on drift network
    "velocity_guidance_weight",  # weight on velocity soft-constraint (0 = off)
    "smoothness_weight",         # regulariser on ||∂_x v||
    "score_matching_weight",     # weight on score-matching objective
)


def prior_sensitivity_report(
    results_by_prior: dict[str, dict],
    knobs: tuple[str, ...] = PRIOR_KNOBS,
    ground_truth_peak: Optional[float] = None,
    peak_window: float = 0.1,
) -> str:
    """
    Sweep across PRIOR configurations (not model class) and report the
    envelope of each output quantity. This is the identifiability probe:
    conclusions that swing under reasonable prior variation are not
    identifiable from snapshot density alone; they are prior-selected.

    Parameters
    ----------
    results_by_prior : dict[str, dict]
        Mapping from prior-config label to a dict containing:
            - "prior"              : dict of knob → value used for this run
            - "archetypes"         : Tensor (K, D, D)              (optional)
            - "instability_curve"  : np.ndarray (T,)               (optional)
            - "regime_fractions"   : dict[str, float]              (optional)
            - "auroc"              : float                         (optional)
        At least two configs required. Configs that differ only in a
        non-prior knob (depth, hidden, seed) will be flagged: they do not
        contribute to identifiability accounting.
    knobs : tuple[str, ...]
        Which knobs are considered part of the prior for this report.
    ground_truth_peak, peak_window
        Forwarded to ``instability_peak_overlap`` to score accuracy of
        peak timing against a known reference.

    Returns
    -------
    str
        Formatted report with per-quantity precision/accuracy envelopes
        and an explicit identifiable/non-identifiable verdict per quantity
        (identifiable = within-tolerance across the prior family).
    """
    configs = list(results_by_prior.keys())
    if len(configs) < 2:
        return (
            "prior_sensitivity_report: need >= 2 prior configurations to probe "
            "identifiability. Vary at least one of: " + ", ".join(knobs)
        )

    lines = [
        "Prior sensitivity report (identifiability probe)",
        "=" * 60,
        f"Prior knobs varied over: {', '.join(knobs)}",
        "",
    ]

    # Flag any config whose 'prior' dict is missing declared knobs — a common
    # mistake is to pass architecture sweeps to this function.
    all_priors = [results_by_prior[c].get("prior", {}) for c in configs]
    varied = {k for k in knobs if len({p.get(k) for p in all_priors}) > 1}
    fixed = {k for k in knobs if k not in varied and any(k in p for p in all_priors)}
    if not varied:
        lines.append(
            "  WARNING: no listed prior knob varies across configs. This "
            "report will not measure identifiability."
        )
    else:
        lines.append(f"  Actually varied: {sorted(varied)}")
    if fixed:
        lines.append(f"  Held fixed:      {sorted(fixed)}")
    lines.append("=" * 60)

    # --- Archetypes: cosine similarity across the prior family -----------
    if all("archetypes" in results_by_prior[c] for c in configs):
        ref = results_by_prior[configs[0]]["archetypes"]
        sims = []
        lines.append("\nArchetype cosine similarity vs first prior:")
        for label in configs[1:]:
            sim = archetype_cosine_similarity(ref, results_by_prior[label]["archetypes"])
            sims.append(sim["median"])
            lines.append(
                f"  {label}: median = {sim['median']:.3f}, min = {sim['min']:.3f}"
            )
        min_sim = min(sims) if sims else 1.0
        verdict = (
            "IDENTIFIABLE (given prior family)"
            if min_sim >= 0.90
            else "NOT IDENTIFIABLE — swings under prior variation"
        )
        lines.append(f"  Verdict (archetype direction): {verdict}  [min={min_sim:.3f}]")

    # --- Instability peak timing: precision AND accuracy -----------------
    if all("instability_curve" in results_by_prior[c] for c in configs):
        curves = [results_by_prior[c]["instability_curve"] for c in configs]
        pk = instability_peak_overlap(
            curves, peak_window=peak_window, ground_truth_peak=ground_truth_peak
        )
        lines.append("\nInstability peak timing:")
        for label, loc in zip(configs, pk["peak_locations"]):
            lines.append(f"  {label}: peak at t = {loc:.3f}")
        lines.append(f"  Std across priors: {pk['peak_std']:.4f}")
        lines.append(
            f"  Cross-prior precision (pairs within {peak_window}): "
            f"{pk['overlap_frac']:.3f}"
        )
        if ground_truth_peak is not None:
            lines.append(
                f"  Ground-truth accuracy: mean |Δt| = {pk['gt_mean_error']:.4f}, "
                f"frac within {peak_window} = {pk['gt_accuracy']:.3f}"
            )
            verdict = (
                "IDENTIFIABLE (accurate under prior variation)"
                if pk["gt_accuracy"] >= 0.80 and pk["peak_std"] <= peak_window
                else "NOT IDENTIFIABLE — magnitude/timing is prior-selected"
            )
        else:
            verdict = (
                "PRECISE but accuracy unknown (pass ground_truth_peak to score)"
                if pk["peak_std"] <= peak_window
                else "NOT IDENTIFIABLE — timing swings under prior variation"
            )
        lines.append(f"  Verdict (peak timing): {verdict}")

    # --- Regime fractions: are 'unstable' / 'plastic' calls prior-stable? -
    if all("regime_fractions" in results_by_prior[c] for c in configs):
        regime_keys = sorted(
            {k for c in configs for k in results_by_prior[c]["regime_fractions"]}
        )
        lines.append("\nRegime fractions across prior family:")
        max_swing = 0.0
        for rk in regime_keys:
            vals = [
                float(results_by_prior[c]["regime_fractions"].get(rk, 0.0))
                for c in configs
            ]
            swing = max(vals) - min(vals)
            max_swing = max(max_swing, swing)
            lines.append(
                f"  {rk:14s}: range [{min(vals):.3f}, {max(vals):.3f}]  Δ={swing:.3f}"
            )
        verdict = (
            "IDENTIFIABLE"
            if max_swing <= 0.10
            else "NOT IDENTIFIABLE — regime calls swing >0.10 under prior variation"
        )
        lines.append(f"  Verdict (regime fractions): {verdict}  [max Δ={max_swing:.3f}]")

    # --- AUROC ------------------------------------------------------------
    aurocs = {
        c: results_by_prior[c]["auroc"]
        for c in configs
        if "auroc" in results_by_prior[c]
    }
    if aurocs:
        lines.append("\nAUROC across priors:")
        for label, auroc in aurocs.items():
            lines.append(f"  {label}: {auroc:.4f}")
        vals = list(aurocs.values())
        delta = max(vals) - min(vals)
        lines.append(f"  Range: {min(vals):.4f} – {max(vals):.4f}  (ΔAUROC = {delta:.4f})")

    lines.append("")
    lines.append(
        "Note: identifiability verdicts are conditional on the declared prior "
        "family. Quantities marked NOT IDENTIFIABLE require an external "
        "constraint (metabolic labelling, timed snapshots, lineage barcoding) "
        "before absolute interpretation."
    )
    return "\n".join(lines)
