"""
Operator Regime Classifier

Classifies cells into operator regimes based on eigenvalue-derived metrics.
"""

import numpy as np
import warnings
from typing import Dict, Tuple, Optional


# Legacy hand-set thresholds. Kept for reproducibility of prior runs — pass
# thresholds='legacy' to opt in. They collapse the 'plastic' class at high
# dim (see class docstring / Task 2.2 notes).
_LEGACY_THRESHOLDS = dict(
    threshold_unstable=0.1,
    threshold_plastic=0.05,
    threshold_deeply_stable=-1.0,
    plasticity_threshold=0.3,
)


class OperatorRegimeClassifier:
    """
    Classify cells into operator regimes based on stability metrics.

    Operator regimes (mutually exclusive by construction):

    - **Unstable**: λ_max⁺ > τ_upper — transition states, near-bifurcation.
    - **Plastic**:  λ_max⁺ in the near-zero band — progenitor states,
      decision points. Defined by BAND membership rather than
      precedence-after-not-unstable-and-not-deeply-stable, so it can no
      longer be silently collapsed by a too-aggressive deep-stable rule.
    - **Deeply Stable**: λ_max⁺ ≤ τ_upper AND stability_depth_q ≤
      deep_stable_cutoff. The cutoff is normalized by the dataset's own
      spectral scale (median |Re(λ)|) so it stays comparable across dim.
    - **Stable**: everything else (λ_max⁺ ≤ τ_upper AND not deeply stable).

    Two structural fixes vs the legacy classifier:

    (a) `stability_depth_q` (per-cell 5th-percentile of Re(λ), computed by
        OperatorMetrics) replaces min(Re(λ)). Extrema shrink with dim purely
        from having more eigenvalues to sample; quantiles don't. Under the
        legacy code at dim=50 nearly every non-unstable cell landed in
        deeply-stable, which is why the plastic class collapsed.

    (b) Thresholds default to data-driven. The unstable / plastic /
        stable-side split comes from terciles of λ_max⁺; the deep-stable
        cutoff comes from k * median(|Re(λ)|). Fully reproducible on the
        same input; comparable across dim. Pass explicit numeric values or
        thresholds='legacy' to override.

    Args:
        thresholds: 'auto' (default, data-driven), 'legacy' (old hardcoded
            0.1 / 0.05 / -1.0 / 0.3), or a dict of override values whose
            keys match the legacy kwargs.
        threshold_unstable: Absolute override for τ_upper. If given, wins
            over 'auto'.
        threshold_plastic: Absolute override for the half-width of the
            near-zero (plastic) band around λ_max⁺ = 0. If given, wins
            over 'auto'.
        threshold_deeply_stable: Absolute override for the deep-stable
            cutoff on stability_depth_q. If given, wins over 'auto'.
        plasticity_threshold: Minimum plasticity index for the plastic
            band. Cells inside the near-zero λ_max⁺ band with plasticity
            below this land in Stable, not Plastic. Default 0.3 (kept —
            plasticity is already dimension-normalized as a fraction).
        deep_stable_scale_k: When thresholds='auto', deep_stable_cutoff =
            -deep_stable_scale_k * median(|Re(λ)|). Default 2.0.
        unstable_tercile: When thresholds='auto', τ_upper = quantile of
            λ_max⁺ at this level. Default 2/3.
    """

    def __init__(
        self,
        thresholds: str = 'hyperbolic',
        threshold_unstable: Optional[float] = None,
        threshold_plastic: Optional[float] = None,
        threshold_deeply_stable: Optional[float] = None,
        plasticity_threshold: float = 0.3,
        deep_stable_scale_k: float = 2.0,
        plastic_scale_k: float = 0.3,
        unstable_tercile: float = 2.0 / 3.0,
    ):
        # thresholds mode:
        #   'hyperbolic' (default, recommended): absolute τ_upper = 0. A cell
        #       is unstable iff its Jacobian has any strictly positive real
        #       eigenvalue. Aligns with the gauge-invariance argument in the
        #       docstring — sign(Re(λ_max)) is the topologically-robust
        #       quantity; the secondary axes (deep-stable cutoff, near-zero
        #       band) are scaled by the dataset's own spectral median. Cross-
        #       dataset and cross-cohort comparisons are meaningful.
        #   'rank': historical tercile split — τ_upper = quantile_{2/3}(λ_max⁺)
        #       so exactly 33% of cells land in unstable BY CONSTRUCTION.
        #       Kept for reproducing prior runs and for use cases that
        #       explicitly want a within-dataset ranking, NOT the default.
        #       'auto' is retained as an alias for 'rank' (deprecated
        #       spelling; will be removed).
        #   'legacy': old hard-coded 0.1 / 0.05 / -1.0 / 0.3 constants.
        if isinstance(thresholds, dict):
            base = dict(_LEGACY_THRESHOLDS); base.update(thresholds)
            self._mode = 'manual'
        elif thresholds == 'legacy':
            base = dict(_LEGACY_THRESHOLDS)
            self._mode = 'legacy'
        elif thresholds == 'hyperbolic':
            base = dict(threshold_unstable=None, threshold_plastic=None,
                        threshold_deeply_stable=None,
                        plasticity_threshold=plasticity_threshold)
            self._mode = 'hyperbolic'
        elif thresholds in ('rank', 'auto'):
            base = dict(threshold_unstable=None, threshold_plastic=None,
                        threshold_deeply_stable=None,
                        plasticity_threshold=plasticity_threshold)
            self._mode = 'rank'
        else:
            raise ValueError(
                f"thresholds must be 'hyperbolic', 'rank', 'legacy', or a dict; "
                f"got {thresholds!r}"
            )

        # Explicit numeric kwargs win regardless of mode.
        if threshold_unstable is not None:
            base['threshold_unstable'] = threshold_unstable
        if threshold_plastic is not None:
            base['threshold_plastic'] = threshold_plastic
        if threshold_deeply_stable is not None:
            base['threshold_deeply_stable'] = threshold_deeply_stable
        base['plasticity_threshold'] = plasticity_threshold

        self._config = base
        self.deep_stable_scale_k = deep_stable_scale_k
        self.plastic_scale_k = plastic_scale_k
        self.unstable_tercile = unstable_tercile

        # Resolved thresholds populated by classify(). classify_with_confidence
        # reuses them so both paths reason about the same numbers.
        self.tau: Optional[float] = None
        self.eps: Optional[float] = None
        self.deeply_stable_threshold: Optional[float] = None
        self.plasticity_threshold: float = plasticity_threshold
        self.spectral_scale: Optional[float] = None

    def _resolve_thresholds(self, metrics: Dict[str, np.ndarray]) -> Dict[str, float]:
        """Fill any None threshold from the metrics distribution.

        Defaults depend on self._mode:

        - hyperbolic (default):
            tau_upper = 0                    (ABSOLUTE — sign of Re(λ_max))
            eps       = plastic_scale_k * spectral_scale
            deep_cut  = -deep_stable_scale_k * spectral_scale

        - rank (legacy behavior, formerly 'auto'):
            tau_upper = quantile_{2/3}(λ_max⁺)   (tercile — pins unstable at 33%)
            eps       = tau_upper / 2
            deep_cut  = -deep_stable_scale_k * spectral_scale

        The spectral scale is median|Re(λ)| across all cells (approximated
        from lambda_max_plus and lambda_min_minus). Explicit numeric
        overrides in self._config always win.
        """
        cfg = dict(self._config)

        lam_max = np.asarray(metrics['lambda_max_plus'])
        lam_min = np.asarray(metrics['lambda_min_minus'])
        # Spectral scale — median absolute Re(λ) across all cells, using
        # both extrema as a proxy for the spectrum's dynamic range. Cheap
        # and does not require store_eigenvalues.
        scale = float(np.median(np.abs(np.concatenate([lam_max, lam_min]))))
        # Guard: constant-spectrum inputs (all zeros) — fall back to a
        # small positive scale so bands don't collapse to zero-width.
        if not np.isfinite(scale) or scale <= 0:
            scale = 1.0
        self.spectral_scale = scale

        if cfg['threshold_unstable'] is None:
            if self._mode == 'hyperbolic':
                cfg['threshold_unstable'] = 0.0
            else:  # rank (or 'auto' alias)
                cfg['threshold_unstable'] = float(np.quantile(lam_max, self.unstable_tercile))

        if cfg['threshold_plastic'] is None:
            if self._mode == 'hyperbolic':
                # Near-zero band width proportional to spectral scale, NOT
                # to tau_upper (which is 0 in hyperbolic mode).
                cfg['threshold_plastic'] = self.plastic_scale_k * scale
            else:  # rank mode: proportional to tau_upper
                cfg['threshold_plastic'] = 0.5 * max(cfg['threshold_unstable'], 1e-6)

        if cfg['threshold_deeply_stable'] is None:
            cfg['threshold_deeply_stable'] = -self.deep_stable_scale_k * scale

        self.tau = cfg['threshold_unstable']
        self.eps = cfg['threshold_plastic']
        self.deeply_stable_threshold = cfg['threshold_deeply_stable']
        return cfg

    def _stability_axis(self, metrics: Dict[str, np.ndarray]) -> np.ndarray:
        """Return the per-cell stability score used for the deep-stable rule.

        Prefers 'stability_depth_q' (dimension-stable quantile); falls back
        to 'lambda_min_minus' (extremum) with a warning for backwards compat.
        """
        if 'stability_depth_q' in metrics:
            return np.asarray(metrics['stability_depth_q'])
        warnings.warn(
            "metrics dict lacks 'stability_depth_q'; falling back to the "
            "dimension-unstable extremum 'lambda_min_minus'. Recompute via "
            "OperatorMetrics from the current codebase to get the quantile-"
            "based stability depth.",
            RuntimeWarning,
            stacklevel=3,
        )
        return np.asarray(metrics['lambda_min_minus'])

    def classify(
        self,
        metrics: Dict[str, np.ndarray]
    ) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
        """
        Classify cells into operator regimes.

        The four regimes are mutually exclusive by construction (built from
        two orthogonal axes — λ_max⁺ band and stability depth — rather than
        the old precedence chain that could silently collapse 'plastic').

        Args:
            metrics: Dictionary with keys 'lambda_max_plus', 'lambda_min_minus',
                'stability_depth_q' (preferred; falls back to 'lambda_min_minus'
                with a warning if absent), 'plasticity', 'stable_dim'.

        Returns:
            regimes: Array of regime labels (n_cells,)
            regime_masks: Dictionary of boolean masks for each regime.
                Sum across regimes equals 1 for every cell (partition).
        """
        cfg = self._resolve_thresholds(metrics)
        tau_upper = cfg['threshold_unstable']
        eps = cfg['threshold_plastic']
        deep_cut = cfg['threshold_deeply_stable']
        p_min = cfg['plasticity_threshold']

        lambda_max = np.asarray(metrics['lambda_max_plus'])
        plasticity = np.asarray(metrics['plasticity'])
        stability = self._stability_axis(metrics)

        n_cells = len(lambda_max)
        regimes = np.empty(n_cells, dtype=object)

        # Three-way partition on λ_max⁺, orthogonal to the deep-stable check.
        unstable_mask = lambda_max > tau_upper
        near_zero = (~unstable_mask) & (np.abs(lambda_max) < eps)
        low_lambda = (~unstable_mask) & (~near_zero)  # everything below the plastic band

        # Plastic requires the near-zero band AND enough near-neutral modes;
        # otherwise it falls through to Stable.
        plastic_mask = near_zero & (plasticity > p_min)
        near_zero_but_not_plastic = near_zero & ~plastic_mask

        # Deep-stable only competes for the low-lambda + near-zero-but-not-
        # plastic slice — never contests the unstable or genuinely-plastic
        # slices.
        candidate_stable = low_lambda | near_zero_but_not_plastic
        deeply_stable_mask = candidate_stable & (stability < deep_cut)
        stable_mask = candidate_stable & ~deeply_stable_mask

        regimes[unstable_mask] = 'unstable'
        regimes[plastic_mask] = 'plastic'
        regimes[deeply_stable_mask] = 'deeply_stable'
        regimes[stable_mask] = 'stable'

        regime_masks = {
            'unstable': unstable_mask,
            'plastic': plastic_mask,
            'deeply_stable': deeply_stable_mask,
            'stable': stable_mask,
        }

        # Partition invariant: exactly one True per cell.
        stacked_sum = sum(m.astype(np.int8) for m in regime_masks.values())
        if not np.all(stacked_sum == 1):
            raise AssertionError(
                "regime_masks do not partition the cells (this is a bug in "
                "the classifier, not a data issue)."
            )

        return regimes, regime_masks

    def sweep_thresholds(
        self,
        metrics: Dict[str, np.ndarray],
        tau_grid: Optional[np.ndarray] = None,
        deep_grid: Optional[np.ndarray] = None,
        n_grid: int = 11,
    ) -> Dict[str, np.ndarray]:
        """
        Regime-fraction sensitivity to (tau_upper, deep_stable_cutoff).

        Produces the numbers for the threshold-sensitivity supplementary
        figure the work order asks for as pre-emptive review defense.
        Rendering is deliberately left to the caller — the library has no
        matplotlib coupling here.

        Args:
            metrics: Same metrics dict as classify().
            tau_grid: Overrides for τ_upper. Defaults to a linspace across
                the observed λ_max⁺ range.
            deep_grid: Overrides for the deep-stable cutoff. Defaults to a
                linspace across the observed stability_depth_q range.
            n_grid: Grid size when the arrays are not supplied.

        Returns:
            Dict with keys:
                'tau_grid', 'deep_grid': the axis values, length T and D.
                For each regime name ('unstable', 'plastic', 'deeply_stable',
                'stable'): a (D, T) array of fractions.
                'reference': dict of the resolved auto thresholds so the
                    caller can overlay a marker on the figure.
        """
        # Trigger threshold resolution so 'reference' is populated.
        _regimes, _masks = self.classify(metrics)
        reference = dict(
            tau_upper=self.tau,
            deep_cut=self.deeply_stable_threshold,
            eps=self.eps,
            spectral_scale=self.spectral_scale,
        )

        lam_max = np.asarray(metrics['lambda_max_plus'])
        stability = self._stability_axis(metrics)

        if tau_grid is None:
            tau_grid = np.linspace(
                float(np.quantile(lam_max, 0.05)),
                float(np.quantile(lam_max, 0.95)),
                n_grid,
            )
        if deep_grid is None:
            deep_grid = np.linspace(
                float(np.quantile(stability, 0.05)),
                float(np.quantile(stability, 0.95)),
                n_grid,
            )

        base = dict(self._config)
        base['plasticity_threshold'] = self.plasticity_threshold

        regimes_axes = ('unstable', 'plastic', 'deeply_stable', 'stable')
        out = {name: np.zeros((len(deep_grid), len(tau_grid))) for name in regimes_axes}

        for i, dcut in enumerate(deep_grid):
            for j, tau in enumerate(tau_grid):
                cls = OperatorRegimeClassifier(
                    thresholds='auto',
                    threshold_unstable=float(tau),
                    threshold_deeply_stable=float(dcut),
                    plasticity_threshold=self.plasticity_threshold,
                    deep_stable_scale_k=self.deep_stable_scale_k,
                    unstable_tercile=self.unstable_tercile,
                )
                _, m = cls.classify(metrics)
                for name in regimes_axes:
                    out[name][i, j] = m[name].mean()

        out['tau_grid'] = np.asarray(tau_grid)
        out['deep_grid'] = np.asarray(deep_grid)
        out['reference'] = reference
        return out
    
    def classify_with_confidence(
        self,
        metrics: Dict[str, np.ndarray]
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Classify cells and return a per-cell BOUNDARY DISTANCE — not a
        probability.

        The returned score is `distance from the decision boundary /
        boundary threshold`, clipped to [0, 1]. It is not calibrated, it
        has no probabilistic interpretation, and it scales inversely with
        an arbitrary threshold constant. Reading it as a posterior
        probability is wrong.

        Task 2.3 renamed the corresponding adata.obs column from
        `regime_confidence` to `boundary_distance` for this reason. The
        method name is preserved here for backwards-compatible callsites;
        callers should treat the returned array as a *distance*, not a
        *confidence*, until Task 3.2 introduces ensemble-concordance
        confidence (which IS a real uncertainty quantifier).

        Args:
            metrics: Dictionary of operator metrics

        Returns:
            regimes: Array of regime labels (n_cells,)
            boundary_distance: Per-cell scores (n_cells,) in [0, 1]. Named
                `confidence` in the local variable and prior return-name for
                back-compat; do not interpret as a probability.
        """
        # classify() resolves and stores self.tau / self.deeply_stable_threshold
        # so the confidence path reasons about the same numbers.
        regimes, _ = self.classify(metrics)

        lambda_max = np.asarray(metrics['lambda_max_plus'])
        stability = self._stability_axis(metrics)
        plasticity = np.asarray(metrics['plasticity'])

        # Guard against divide-by-zero if a user picks tau=0.
        tau_denom = max(abs(self.tau), 1e-6)
        deep_denom = max(abs(self.deeply_stable_threshold), 1e-6)

        n_cells = len(regimes)
        confidence = np.zeros(n_cells)

        for i in range(n_cells):
            if regimes[i] == 'unstable':
                confidence[i] = np.clip((lambda_max[i] - self.tau) / tau_denom, 0, 1)
            elif regimes[i] == 'deeply_stable':
                confidence[i] = np.clip(
                    (self.deeply_stable_threshold - stability[i]) / deep_denom, 0, 1
                )
            elif regimes[i] == 'plastic':
                confidence[i] = np.clip(
                    (plasticity[i] - self.plasticity_threshold)
                    / max(1 - self.plasticity_threshold, 1e-6),
                    0, 1,
                )
            else:  # stable
                confidence[i] = np.clip((self.tau - lambda_max[i]) / tau_denom, 0, 1)

        return regimes, confidence
    
    def get_regime_statistics(
        self,
        regimes: np.ndarray,
        metrics: Dict[str, np.ndarray]
    ) -> Dict[str, Dict[str, float]]:
        """
        Compute summary statistics for each regime.
        
        Args:
            regimes: Array of regime labels
            metrics: Dictionary of operator metrics
            
        Returns:
            Dictionary mapping regime names to statistics
        """
        unique_regimes = np.unique(regimes)
        stats = {}
        
        for regime in unique_regimes:
            mask = regimes == regime
            
            stats[regime] = {
                'count': mask.sum(),
                'fraction': mask.mean(),
                'lambda_max_mean': metrics['lambda_max_plus'][mask].mean(),
                'lambda_max_std': metrics['lambda_max_plus'][mask].std(),
                'lambda_min_mean': metrics['lambda_min_minus'][mask].mean(),
                'lambda_min_std': metrics['lambda_min_minus'][mask].std(),
                'plasticity_mean': metrics['plasticity'][mask].mean(),
                'plasticity_std': metrics['plasticity'][mask].std(),
                'stable_dim_mean': metrics['stable_dim'][mask].mean(),
                'stable_dim_std': metrics['stable_dim'][mask].std(),
            }
        
        return stats
    
    def compare_regimes_across_conditions(
        self,
        regimes: np.ndarray,
        conditions: np.ndarray,
        celltype: Optional[np.ndarray] = None
    ) -> Dict[str, Dict]:
        """
        Compare operator regime distributions across conditions.
        
        This is critical for showing non-redundancy with expression-based cell types.
        
        Args:
            regimes: Array of regime labels (n_cells,)
            conditions: Array of condition labels (n_cells,)
            celltype: Optional array of cell type labels (n_cells,)
            
        Returns:
            Dictionary with comparison statistics
        """
        unique_conditions = np.unique(conditions)
        comparison = {}
        
        # Overall regime distribution per condition
        for cond in unique_conditions:
            cond_mask = conditions == cond
            cond_regimes = regimes[cond_mask]
            
            regime_counts = {}
            for regime in ['stable', 'plastic', 'unstable', 'deeply_stable']:
                regime_counts[regime] = (cond_regimes == regime).sum()
            
            comparison[cond] = {
                'total_cells': cond_mask.sum(),
                'regime_counts': regime_counts,
                'regime_fractions': {
                    k: v / cond_mask.sum() for k, v in regime_counts.items()
                }
            }
        
        # If cell types provided, compare within cell types
        if celltype is not None:
            comparison['by_celltype'] = {}
            unique_celltypes = np.unique(celltype)
            
            for ct in unique_celltypes:
                ct_mask = celltype == ct
                comparison['by_celltype'][ct] = {}
                
                for cond in unique_conditions:
                    cond_ct_mask = ct_mask & (conditions == cond)
                    if cond_ct_mask.sum() == 0:
                        continue
                    
                    cond_ct_regimes = regimes[cond_ct_mask]
                    
                    regime_counts = {}
                    for regime in ['stable', 'plastic', 'unstable', 'deeply_stable']:
                        regime_counts[regime] = (cond_ct_regimes == regime).sum()
                    
                    comparison['by_celltype'][ct][cond] = {
                        'total_cells': cond_ct_mask.sum(),
                        'regime_counts': regime_counts,
                        'regime_fractions': {
                            k: v / cond_ct_mask.sum() if cond_ct_mask.sum() > 0 else 0
                            for k, v in regime_counts.items()
                        }
                    }
        
        return comparison
