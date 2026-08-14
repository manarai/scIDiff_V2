"""
Statistical inference for operator-regime × condition comparisons.

Task 3.1. Before this module, `atlas.validate_nonredundancy` returned regime
counts, per-cell-type regime diversity, and a Shannon entropy — no null
model, no p-value, no effect size, no confidence interval, no multiple-
testing correction. "Naïve T cells are 60% deeply-stable in old and 45% in
young" is not a result without a null distribution, and the plug-in Shannon
estimator it used is downward-biased at small n so uneven sampling depth
alone could manufacture a signal.

This module provides the missing pieces:

- Chi-square + Cramér's V on regime × condition tables per cell type
  (`chi_square_regime_vs_condition`).
- Permutation null that shuffles CONDITION labels WITHIN cell type
  (`permutation_p_value`) — preserves cell-type composition and regime
  marginals while destroying the condition association.
- Depth-matched subsampling that equalizes per-(celltype × condition) both
  cell count AND UMI depth via quantile bins
  (`depth_matched_subsample`), then a distribution of statistics across
  many matched draws (`matched_test`).
- Benjamini–Hochberg FDR across cell types (`bh_correction`).
- Bootstrap CIs on regime fractions (`bootstrap_ci_regime_fractions`).
- Miller-Madow entropy bias correction (`miller_madow_entropy`) for
  callers who still want an entropy summary.
- `test_nonredundancy(...)` — the high-level entry point that ties the
  above together and returns a single structured dict.

Design note. The critical decision the whole task hinges on is depth
matching: regime assignment depends on the drift field, which depends on
the embedding, which depends on UMI depth and batch. If a "condition
effect" vanishes under matching, the effect was depth. Callers who cannot
supply counts get a warning and a matched-count-only run; the returned
dict marks depth_matching_applied so downstream code can see it.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.stats import chi2_contingency


ArrayLike = np.ndarray
REGIME_NAMES: Tuple[str, ...] = ("stable", "plastic", "unstable", "deeply_stable")


# -----------------------------------------------------------------------------
# Effect size and single-sample statistic
# -----------------------------------------------------------------------------

def _contingency_table(
    regimes: ArrayLike,
    conditions: ArrayLike,
    regime_names: Sequence[str] = REGIME_NAMES,
) -> Tuple[np.ndarray, List[str]]:
    """Return a (n_regimes, n_conditions) integer table.

    Zeros are kept in place so shapes are stable across cell types.
    """
    conds = list(np.unique(conditions))
    table = np.zeros((len(regime_names), len(conds)), dtype=np.int64)
    for i, r in enumerate(regime_names):
        for j, c in enumerate(conds):
            table[i, j] = int(((regimes == r) & (conditions == c)).sum())
    return table, conds


def chi_square_regime_vs_condition(
    regimes: ArrayLike,
    conditions: ArrayLike,
    regime_names: Sequence[str] = REGIME_NAMES,
) -> Dict[str, float]:
    """Chi-square + Cramér's V on a regime × condition contingency table.

    Cramér's V is a [0, 1] effect size that's directly comparable across
    cell types with different N — you must report it alongside the p-value
    because chi-square p-values become nominally significant at any
    trivial effect once N grows.

    Returns:
        dict with keys chi2, p, dof, cramers_v, n_cells, table (2-D int),
        conditions (list).
    """
    table, conds = _contingency_table(regimes, conditions, regime_names)
    n = int(table.sum())
    # Drop all-zero rows so chi2_contingency doesn't warn/fail; the effect
    # size denominator uses ORIGINAL n (all cells that were in the input).
    nonzero_rows = table.sum(axis=1) > 0
    reduced = table[nonzero_rows]

    if reduced.shape[0] < 2 or reduced.shape[1] < 2 or n == 0:
        # Undefined chi-square (needs ≥ 2 rows and ≥ 2 cols with data).
        return dict(
            chi2=float("nan"), p=float("nan"), dof=0, cramers_v=float("nan"),
            n_cells=n, table=table, conditions=conds,
        )

    chi2, p, dof, _expected = chi2_contingency(reduced)
    r, c = reduced.shape
    denom = n * (min(r, c) - 1)
    cramers_v = float(np.sqrt(chi2 / denom)) if denom > 0 else float("nan")
    return dict(
        chi2=float(chi2), p=float(p), dof=int(dof), cramers_v=cramers_v,
        n_cells=n, table=table, conditions=conds,
    )


# -----------------------------------------------------------------------------
# Permutation null
# -----------------------------------------------------------------------------

def permutation_p_value(
    regimes: ArrayLike,
    conditions: ArrayLike,
    n_perm: int = 1000,
    rng: Optional[np.random.Generator] = None,
    regime_names: Sequence[str] = REGIME_NAMES,
) -> Dict[str, float]:
    """Permutation p-value for the regime ↔ condition association.

    Shuffles CONDITION labels (not regime labels), which preserves both
    cell-type composition and the regime marginal while destroying the
    condition-effect signal. Reported p is the standard (1 + count) / (1
    + n_perm) form — never zero, so log-plots don't break.

    Returns:
        dict with observed_chi2, p_perm, cramers_v (observed), null_mean,
        null_std, n_perm_effective.
    """
    if rng is None:
        rng = np.random.default_rng(0)

    observed = chi_square_regime_vs_condition(regimes, conditions, regime_names)
    obs_chi2 = observed["chi2"]
    if not np.isfinite(obs_chi2):
        return dict(
            observed_chi2=float("nan"), p_perm=float("nan"),
            cramers_v=float("nan"), null_mean=float("nan"),
            null_std=float("nan"), n_perm_effective=0,
        )

    regimes = np.asarray(regimes)
    conditions = np.asarray(conditions)
    n_greater = 0
    null_chi2 = np.empty(n_perm, dtype=float)
    n_used = 0
    for k in range(n_perm):
        shuffled = rng.permutation(conditions)
        stat = chi_square_regime_vs_condition(regimes, shuffled, regime_names)["chi2"]
        if not np.isfinite(stat):
            continue
        null_chi2[n_used] = stat
        n_used += 1
        if stat >= obs_chi2:
            n_greater += 1

    if n_used == 0:
        return dict(
            observed_chi2=obs_chi2, p_perm=float("nan"),
            cramers_v=observed["cramers_v"], null_mean=float("nan"),
            null_std=float("nan"), n_perm_effective=0,
        )

    p_perm = (1 + n_greater) / (1 + n_used)
    null_chi2 = null_chi2[:n_used]
    return dict(
        observed_chi2=obs_chi2,
        p_perm=float(p_perm),
        cramers_v=observed["cramers_v"],
        null_mean=float(null_chi2.mean()),
        null_std=float(null_chi2.std(ddof=1)) if n_used > 1 else 0.0,
        n_perm_effective=int(n_used),
    )


# -----------------------------------------------------------------------------
# Depth-matched subsampling
# -----------------------------------------------------------------------------

def depth_matched_subsample(
    celltypes: ArrayLike,
    conditions: ArrayLike,
    counts: Optional[ArrayLike],
    n_bins: int = 5,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    """Return a boolean mask that equalizes (celltype × condition) N and
    UMI-depth distribution.

    For each cell type:
      1. Bin cells by UMI-depth quantiles (`n_bins` equal-frequency bins).
      2. Within each bin, downsample every condition to the minimum count
         so all conditions contribute equally to that bin.
      3. Retain only the intersected mask.

    If `counts` is None, degrade to N-only matching (equal cells per
    condition per celltype) and issue a RuntimeWarning — the caller then
    cannot rule out depth as the confounder.

    Returns:
        Boolean mask of length n_cells.
    """
    if rng is None:
        rng = np.random.default_rng(0)

    celltypes = np.asarray(celltypes)
    conditions = np.asarray(conditions)
    n = len(celltypes)
    keep = np.zeros(n, dtype=bool)

    if counts is None:
        warnings.warn(
            "depth_matched_subsample called without counts — falling back "
            "to N-only matching. Depth remains a potential confounder; the "
            "caller cannot rule out that observed effects are depth "
            "artifacts.",
            RuntimeWarning,
            stacklevel=2,
        )
    else:
        counts = np.asarray(counts, dtype=float)

    for ct in np.unique(celltypes):
        ct_mask = celltypes == ct
        ct_idx = np.flatnonzero(ct_mask)
        ct_conds = conditions[ct_mask]
        uniq_conds = np.unique(ct_conds)

        if len(uniq_conds) < 2:
            # Nothing to match against; keep all cells of this celltype.
            keep[ct_idx] = True
            continue

        if counts is None:
            # N-only match: subsample every condition to min N.
            per_cond = [ct_idx[ct_conds == c] for c in uniq_conds]
            n_keep = min(len(x) for x in per_cond)
            for arr in per_cond:
                sel = rng.choice(arr, size=n_keep, replace=False)
                keep[sel] = True
            continue

        # Depth-matched. Quantile bins on this celltype's counts.
        ct_counts = counts[ct_mask]
        # Use unique-quantile edges so extreme ties don't break qcut.
        edges = np.quantile(ct_counts, np.linspace(0, 1, n_bins + 1))
        edges = np.unique(edges)
        if len(edges) < 2:
            # Constant counts → equivalent to N-only within this celltype.
            per_cond = [ct_idx[ct_conds == c] for c in uniq_conds]
            n_keep = min(len(x) for x in per_cond)
            for arr in per_cond:
                sel = rng.choice(arr, size=n_keep, replace=False)
                keep[sel] = True
            continue

        # Bin index per cell (0..len(edges)-2).
        bin_idx = np.clip(
            np.searchsorted(edges, ct_counts, side="right") - 1,
            0, len(edges) - 2,
        )

        for b in range(len(edges) - 1):
            bin_cells = ct_idx[bin_idx == b]
            if bin_cells.size == 0:
                continue
            bin_conds = conditions[bin_cells]
            per_cond = [bin_cells[bin_conds == c] for c in uniq_conds]
            n_keep = min(len(x) for x in per_cond)
            if n_keep == 0:
                continue
            for arr in per_cond:
                sel = rng.choice(arr, size=n_keep, replace=False)
                keep[sel] = True

    return keep


def matched_test(
    regimes: ArrayLike,
    conditions: ArrayLike,
    celltypes: ArrayLike,
    counts: Optional[ArrayLike],
    n_matched: int = 100,
    n_perm: int = 200,
    n_bins: int = 5,
    rng: Optional[np.random.Generator] = None,
    regime_names: Sequence[str] = REGIME_NAMES,
) -> Dict[str, Dict]:
    """Per-cell-type matched permutation test.

    For each of `n_matched` random depth-matched subsamples, runs the
    permutation test per cell type. Reports the DISTRIBUTION of statistics
    and p-values, not a single draw — a single lucky draw is not evidence.

    Returns:
        {celltype: {
            p_perm_median, p_perm_lo, p_perm_hi,
            cramers_v_median, cramers_v_lo, cramers_v_hi,
            n_matched_effective, n_perm_per_draw,
        }, ...}
        plus 'depth_matching_applied': bool.
    """
    if rng is None:
        rng = np.random.default_rng(0)

    regimes = np.asarray(regimes)
    conditions = np.asarray(conditions)
    celltypes = np.asarray(celltypes)

    unique_cts = list(np.unique(celltypes))
    p_draws: Dict[str, List[float]] = {ct: [] for ct in unique_cts}
    v_draws: Dict[str, List[float]] = {ct: [] for ct in unique_cts}

    for k in range(n_matched):
        seed = int(rng.integers(0, 2**31 - 1))
        keep = depth_matched_subsample(
            celltypes, conditions, counts, n_bins=n_bins,
            rng=np.random.default_rng(seed),
        )
        for ct in unique_cts:
            m = keep & (celltypes == ct)
            if m.sum() < 4 or len(np.unique(conditions[m])) < 2:
                continue
            perm_res = permutation_p_value(
                regimes[m], conditions[m], n_perm=n_perm,
                rng=np.random.default_rng(seed + 1),
                regime_names=regime_names,
            )
            if np.isfinite(perm_res["p_perm"]):
                p_draws[ct].append(perm_res["p_perm"])
                v_draws[ct].append(perm_res["cramers_v"])

    out: Dict[str, Dict] = {}
    for ct in unique_cts:
        ps = np.asarray(p_draws[ct])
        vs = np.asarray(v_draws[ct])
        if len(ps) == 0:
            out[ct] = dict(
                p_perm_median=float("nan"), p_perm_lo=float("nan"),
                p_perm_hi=float("nan"),
                cramers_v_median=float("nan"), cramers_v_lo=float("nan"),
                cramers_v_hi=float("nan"),
                n_matched_effective=0, n_perm_per_draw=n_perm,
            )
            continue
        out[ct] = dict(
            p_perm_median=float(np.median(ps)),
            p_perm_lo=float(np.quantile(ps, 0.025)),
            p_perm_hi=float(np.quantile(ps, 0.975)),
            cramers_v_median=float(np.median(vs)),
            cramers_v_lo=float(np.quantile(vs, 0.025)),
            cramers_v_hi=float(np.quantile(vs, 0.975)),
            n_matched_effective=int(len(ps)),
            n_perm_per_draw=n_perm,
        )
    out["_depth_matching_applied"] = counts is not None
    return out


# -----------------------------------------------------------------------------
# Multiple testing
# -----------------------------------------------------------------------------

def bh_correction(p_values: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg FDR-corrected q-values.

    Standard step-up: rank ascending, q_i = min_{k >= i} min(1, p_(k) * m / k).
    NaN inputs pass through as NaN q-values and are excluded from the
    denominator implicitly (their p-values don't rank).
    """
    p = np.asarray(p_values, dtype=float)
    finite = np.isfinite(p)
    if not finite.any():
        return np.full_like(p, np.nan, dtype=float)

    q = np.full_like(p, np.nan, dtype=float)
    p_fin = p[finite]
    m = len(p_fin)
    order = np.argsort(p_fin)
    ranked = p_fin[order]
    # Enforce monotonicity of q from the tail back.
    q_ranked = np.minimum.accumulate((ranked * m / np.arange(1, m + 1))[::-1])[::-1]
    q_ranked = np.clip(q_ranked, 0.0, 1.0)
    q_fin = np.empty_like(q_ranked)
    q_fin[order] = q_ranked
    q[finite] = q_fin
    return q


# -----------------------------------------------------------------------------
# Bootstrap CIs
# -----------------------------------------------------------------------------

def bootstrap_ci_regime_fractions(
    regimes: ArrayLike,
    n_boot: int = 1000,
    ci: float = 0.95,
    rng: Optional[np.random.Generator] = None,
    regime_names: Sequence[str] = REGIME_NAMES,
) -> Dict[str, Tuple[float, float, float]]:
    """Percentile-bootstrap CIs on regime fractions for a single cohort.

    Returns: {regime: (point, lo, hi)}.
    """
    if rng is None:
        rng = np.random.default_rng(0)

    regimes = np.asarray(regimes)
    n = len(regimes)
    point = {r: float((regimes == r).mean()) if n else float("nan") for r in regime_names}
    if n == 0:
        return {r: (float("nan"),) * 3 for r in regime_names}

    frac_boot = np.empty((n_boot, len(regime_names)), dtype=float)
    lo_q, hi_q = (1 - ci) / 2, 1 - (1 - ci) / 2
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        sample = regimes[idx]
        for j, r in enumerate(regime_names):
            frac_boot[b, j] = (sample == r).mean()

    out = {}
    for j, r in enumerate(regime_names):
        lo = float(np.quantile(frac_boot[:, j], lo_q))
        hi = float(np.quantile(frac_boot[:, j], hi_q))
        out[r] = (point[r], lo, hi)
    return out


# -----------------------------------------------------------------------------
# Entropy bias correction
# -----------------------------------------------------------------------------

def miller_madow_entropy(counts: ArrayLike) -> float:
    """Miller-Madow bias-corrected Shannon entropy (bits).

    plug-in H = -sum p log2 p; MM correction adds (K - 1) / (2 n ln 2),
    where K is the number of observed non-zero bins and n is the total
    count. Chao-Shen is more sophisticated but Miller-Madow removes the
    dominant O(1/n) bias term that made the raw entropy estimator inflate
    with sample size — enough for the level of claim scOpAtlas is making.
    """
    counts = np.asarray(counts, dtype=float)
    n = counts.sum()
    if n <= 0:
        return float("nan")
    p = counts / n
    p = p[p > 0]
    plugin = float(-np.sum(p * np.log2(p)))
    K = len(p)
    if K <= 1:
        return plugin
    return plugin + (K - 1) / (2 * n * np.log(2))


# -----------------------------------------------------------------------------
# High-level entry
# -----------------------------------------------------------------------------

@dataclass
class NonredundancyResult:
    """Structured output of test_nonredundancy — attribute access and
    also indexable like a dict (see .to_dict())."""

    per_celltype: Dict[str, Dict] = field(default_factory=dict)
    q_values: Dict[str, float] = field(default_factory=dict)
    matched: Dict[str, Dict] = field(default_factory=dict)
    bootstrap_ci: Dict[str, Dict[str, Dict[str, Tuple[float, float, float]]]] = field(
        default_factory=dict
    )
    depth_matching_applied: bool = False
    n_perm: int = 0
    n_matched: int = 0
    n_boot: int = 0

    def to_dict(self) -> Dict:
        return dict(
            per_celltype=self.per_celltype,
            q_values=self.q_values,
            matched=self.matched,
            bootstrap_ci=self.bootstrap_ci,
            depth_matching_applied=self.depth_matching_applied,
            n_perm=self.n_perm,
            n_matched=self.n_matched,
            n_boot=self.n_boot,
        )


def test_nonredundancy(
    regimes: ArrayLike,
    celltype: ArrayLike,
    condition: ArrayLike,
    counts: Optional[ArrayLike] = None,
    n_perm: int = 1000,
    n_matched: int = 100,
    n_perm_matched: int = 200,
    n_boot: int = 1000,
    n_bins: int = 5,
    rng: Optional[np.random.Generator] = None,
    regime_names: Sequence[str] = REGIME_NAMES,
) -> NonredundancyResult:
    """High-level entry.

    For each cell type, runs:
      1. Per-cell-type chi-square + Cramér's V observed on the full data.
      2. Per-cell-type permutation p-value on the full data.
      3. Per-condition bootstrap CIs on regime fractions.
    Then across cell types:
      4. BH FDR correction on the per-cell-type p-values (q-values).
      5. Matched-permutation test using depth-matched subsamples — the
         acceptance-critical step that separates biology from depth.

    Returns a NonredundancyResult (see its docstring / to_dict()).
    """
    if rng is None:
        rng = np.random.default_rng(0)

    regimes = np.asarray(regimes)
    celltype = np.asarray(celltype)
    condition = np.asarray(condition)
    counts_arr = None if counts is None else np.asarray(counts, dtype=float)

    per_ct: Dict[str, Dict] = {}
    boot_by_ct: Dict[str, Dict[str, Dict[str, Tuple[float, float, float]]]] = {}
    unique_cts = list(np.unique(celltype))

    ordered_ps = []
    for ct in unique_cts:
        m = celltype == ct
        r_ct = regimes[m]
        c_ct = condition[m]

        obs = chi_square_regime_vs_condition(r_ct, c_ct, regime_names)
        perm = permutation_p_value(
            r_ct, c_ct, n_perm=n_perm,
            rng=np.random.default_rng(rng.integers(0, 2**31 - 1)),
            regime_names=regime_names,
        )
        per_ct[ct] = dict(
            n_cells=obs["n_cells"],
            chi2=obs["chi2"],
            p_chi2=obs["p"],
            cramers_v=obs["cramers_v"],
            p_perm=perm["p_perm"],
            null_mean_chi2=perm["null_mean"],
            null_std_chi2=perm["null_std"],
            entropy_mm=miller_madow_entropy(
                _contingency_table(r_ct, c_ct, regime_names)[0].sum(axis=1)
            ),
        )
        ordered_ps.append(perm["p_perm"])

        # Bootstrap CIs per condition within this cell type.
        boot_by_ct[ct] = {}
        for cond in np.unique(c_ct):
            m_cc = m & (condition == cond)
            boot_by_ct[ct][str(cond)] = bootstrap_ci_regime_fractions(
                regimes[m_cc], n_boot=n_boot,
                rng=np.random.default_rng(rng.integers(0, 2**31 - 1)),
                regime_names=regime_names,
            )

    q_arr = bh_correction(ordered_ps)
    q_values = {ct: float(q_arr[i]) for i, ct in enumerate(unique_cts)}

    matched = matched_test(
        regimes, condition, celltype, counts_arr,
        n_matched=n_matched, n_perm=n_perm_matched, n_bins=n_bins,
        rng=np.random.default_rng(rng.integers(0, 2**31 - 1)),
        regime_names=regime_names,
    )
    depth_applied = bool(matched.pop("_depth_matching_applied", False))

    return NonredundancyResult(
        per_celltype=per_ct,
        q_values=q_values,
        matched=matched,
        bootstrap_ci=boot_by_ct,
        depth_matching_applied=depth_applied,
        n_perm=n_perm,
        n_matched=n_matched,
        n_boot=n_boot,
    )


# -----------------------------------------------------------------------------
# Count-source helper for AnnData integration
# -----------------------------------------------------------------------------

def infer_counts_from_adata(adata) -> Tuple[Optional[np.ndarray], str]:
    """Auto-detect per-cell UMI counts on an AnnData.

    Fallback chain (Task 3.1 decision):
      1. adata.obs['n_counts']      (scanpy standard)
      2. adata.obs['total_counts']  (scverse standard)
      3. adata.raw.X row-sums       (if raw is present)
      4. adata.X row-sums           (last resort)

    Returns (counts, source_str). counts is None only if nothing worked.
    """
    for key in ("n_counts", "total_counts"):
        if key in adata.obs.columns:
            return np.asarray(adata.obs[key].values, dtype=float), f"obs['{key}']"

    if getattr(adata, "raw", None) is not None:
        X = adata.raw.X
        try:
            counts = np.asarray(X.sum(axis=1)).ravel().astype(float)
        except Exception:  # pragma: no cover
            counts = None
        if counts is not None and counts.size == adata.n_obs:
            return counts, "raw.X row-sums"

    try:
        counts = np.asarray(adata.X.sum(axis=1)).ravel().astype(float)
    except Exception:  # pragma: no cover
        return None, "unavailable"
    return counts, "X row-sums"
