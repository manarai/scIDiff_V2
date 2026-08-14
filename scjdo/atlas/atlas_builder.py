"""
Stable Operator Atlas Builder

Main interface for building and managing operator-based cell state atlases.
"""

import torch
import numpy as np
import anndata as ad
from typing import Dict, Optional, List, Tuple

from .operator_metrics import OperatorMetrics
from .regime_classifier import OperatorRegimeClassifier
from . import statistics as _stats
from . import quality as _quality


class StableOperatorAtlas:
    """
    Stable Operator Atlas: Defines cellular states by their local stability structure.
    
    This is the main interface for scOpAtlas functionality. It takes a trained
    scJDO drift model and computes operator regimes for all cells in an AnnData object.
    
    Key capabilities:
    - Compute operator metrics (λ_max⁺, λ_min⁻, P, S)
    - Classify cells into operator regimes (stable/plastic/unstable/deeply_stable)
    - Compare regimes across conditions
    - Validate non-redundancy with expression-based cell types
    - Store results in AnnData for downstream analysis
    
    Args:
        adata: AnnData object with single-cell data
        drift_model: Trained DriftField model from scJDO
        use_rep: Key in adata.obsm for state representation (default: "X_pca")
        pseudotime_key: Key in adata.obs for pseudotime (default: "pseudotime")
        device: Device for computation (default: "cpu")
    
    Example:
        >>> from scjdo.models.drift import DriftField
        >>> from scjdo.atlas import StableOperatorAtlas
        >>> 
        >>> # Load data and model
        >>> adata = ad.read_h5ad("data.h5ad")
        >>> drift_model = DriftField.load("my_model.pt")
        >>> 
        >>> # Build atlas
        >>> atlas = StableOperatorAtlas(adata, drift_model)
        >>> atlas.build()
        >>> 
        >>> # Access results
        >>> print(atlas.regimes)
        >>> print(atlas.metrics)
    """
    
    def __init__(
        self,
        adata: ad.AnnData,
        drift_model,
        use_rep: str = "X_pca",
        pseudotime_key: str = "pseudotime",
        device: str = "cpu"
    ):
        self.adata = adata
        self.drift_model = drift_model
        self.use_rep = use_rep
        self.pseudotime_key = pseudotime_key
        self.device = device
        
        # Initialize components
        self.metrics_computer = None
        self.classifier = None
        
        # Results
        self.metrics = None
        self.regimes = None
        self.regime_masks = None
        self.confidence = None
        
        # Validate inputs
        self._validate_inputs()
    
    def _validate_inputs(self):
        """Validate that required data is present."""
        if self.use_rep not in self.adata.obsm:
            raise ValueError(
                f"Representation '{self.use_rep}' not found in adata.obsm. "
                f"Available: {list(self.adata.obsm.keys())}"
            )
        
        if self.pseudotime_key not in self.adata.obs:
            raise ValueError(
                f"Pseudotime key '{self.pseudotime_key}' not found in adata.obs "
                f"(available keys: {list(self.adata.obs.columns)}). "
                "scOpAtlas requires a real pseudotime axis: the drift field's "
                "Jacobian is evaluated at (x, t), and substituting arbitrary "
                "row-order time produces a fully populated but meaningless atlas. "
                "Fix by computing pseudotime first (e.g. sc.tl.dpt) and storing "
                "it in adata.obs, or pass pseudotime_key= pointing at an existing "
                "column."
            )
    
    def build(
        self,
        epsilon: float = 0.1,
        thresholds: str = 'hyperbolic',
        threshold_unstable: Optional[float] = None,
        threshold_plastic: Optional[float] = None,
        threshold_deeply_stable: Optional[float] = None,
        plasticity_threshold: float = 0.3,
        deep_stable_scale_k: float = 2.0,
        batch_size: int = 32,
        compute_confidence: bool = True,
        store_eigenvalues: bool = False,
        quality_check: bool = True,
        force: bool = False,
    ):
        """
        Build the Stable Operator Atlas.

        This is the main method that computes all operator metrics and classifies
        cells into operator regimes.

        Args:
            epsilon: Threshold for near-neutral modes (plasticity index, |λ|)
            thresholds: 'auto' (default) — derive classification thresholds
                from the dataset's own spectrum, so results stay comparable
                across ambient dim. 'legacy' — old hard-coded 0.1 / 0.05 /
                -1.0 / 0.3 (preserved for reproducing prior runs; note the
                dim-collapse of the plastic class documented in the
                classifier).
            threshold_unstable, threshold_plastic, threshold_deeply_stable:
                Explicit numeric overrides. When given, take precedence over
                'thresholds' mode for that specific value.
            plasticity_threshold: Minimum plasticity index for plastic regime.
                Already dimension-normalized as a fraction, so kept fixed.
            batch_size: Batch size for processing
            compute_confidence: Whether to compute confidence scores.
            store_eigenvalues: If True, keep the full per-cell eigenvalue
                spectra and write them to adata.uns['operator_eigenvalues'] as
                complex64. Off by default: at 100k cells × 50 dims the spectrum
                is ~40 MB and travels with every saved .h5ad. Turn on for
                spectrum-diagnostic work.
            quality_check: If True (default), run the drift-quality gate
                before classifying. Raises DriftQualityError on failure;
                the full report always lands in adata.uns['drift_quality'].
                See scjdo.atlas.quality.check_drift_quality.
            force: If True, quality-gate failures print a loud banner and
                classification proceeds anyway. Use only when you've
                understood the failure and have a stated reason.
        """
        print("Building Stable Operator Atlas...")

        if quality_check:
            print("Running drift-quality gate...")
            report = _quality.check_drift_quality(
                model=self.drift_model,
                adata=self.adata,
                use_rep=self.use_rep,
                pseudotime_key=self.pseudotime_key,
                device=self.device,
            )
            self.adata.uns['drift_quality'] = report
            print(report['summary'])
            _quality.raise_if_failed(report, force=force)

        self.metrics_computer = OperatorMetrics(
            drift_model=self.drift_model,
            epsilon=epsilon,
            device=self.device,
        )

        self.classifier = OperatorRegimeClassifier(
            thresholds=thresholds,
            threshold_unstable=threshold_unstable,
            threshold_plastic=threshold_plastic,
            threshold_deeply_stable=threshold_deeply_stable,
            plasticity_threshold=plasticity_threshold,
            deep_stable_scale_k=deep_stable_scale_k,
        )

        X = self.adata.obsm[self.use_rep]
        X_tensor = torch.tensor(X, dtype=torch.float32)

        # Get pseudotime. _validate_inputs already guaranteed the key exists —
        # if it doesn't, construction raised, so we do not fabricate a fallback.
        pseudotime = self.adata.obs[self.pseudotime_key].values

        print("Computing operator metrics...")
        self.metrics = self.metrics_computer.compute_metrics_at_pseudotime(
            X_tensor,
            pseudotime,
            batch_size=batch_size,
            store_eigenvalues=store_eigenvalues,
        )

        print("Classifying operator regimes...")
        # Always populate regime_masks — the confidence path used to leave it
        # None, which broke any downstream code doing atlas.regime_masks[...].
        self.regimes, self.regime_masks = self.classifier.classify(self.metrics)
        if compute_confidence:
            _, self.confidence = self.classifier.classify_with_confidence(self.metrics)
        else:
            self.confidence = None

        self._store_results()

        print("Atlas building complete!")
        self._print_summary()

    def _store_results(self):
        """Store operator metrics and regimes in AnnData object."""
        self.adata.obs['operator_regime'] = self.regimes
        self.adata.obs['lambda_max_plus'] = self.metrics['lambda_max_plus']
        self.adata.obs['lambda_min_minus'] = self.metrics['lambda_min_minus']
        if 'stability_depth_q' in self.metrics:
            self.adata.obs['stability_depth_q'] = self.metrics['stability_depth_q']
        self.adata.obs['plasticity'] = self.metrics['plasticity']
        self.adata.obs['stable_dim'] = self.metrics['stable_dim']

        if self.confidence is not None:
            # Renamed from 'regime_confidence' (Task 2.3): the returned score
            # is a clipped, hand-set-threshold normalized distance to the
            # decision boundary — not a probability, no calibration. Reading
            # it as a posterior probability (which the old name invited) is
            # wrong. The honest per-cell uncertainty quantifier is
            # `regime_concordance` (Task 3.2, OperatorEnsemble), which is
            # the fraction of ensemble members agreeing with the modal
            # regime. Use ensemble concordance when reporting uncertainty;
            # this column is an interim single-model fallback.
            self.adata.obs['boundary_distance'] = self.confidence

        # Only write eigenvalues if the caller opted in via store_eigenvalues.
        if 'eigenvalues' in self.metrics:
            self.adata.uns['operator_eigenvalues'] = self.metrics['eigenvalues']
    
    def _print_summary(self):
        """Print summary statistics of the atlas."""
        print("\n=== Operator Regime Summary ===")
        unique_regimes, counts = np.unique(self.regimes, return_counts=True)
        
        for regime, count in zip(unique_regimes, counts):
            fraction = count / len(self.regimes)
            print(f"{regime:15s}: {count:6d} cells ({fraction*100:5.1f}%)")
        
        print("\n=== Operator Metrics Summary ===")
        for metric_name, metric_values in self.metrics.items():
            if metric_name != 'eigenvalues':
                print(f"{metric_name:20s}: mean={metric_values.mean():7.3f}, std={metric_values.std():7.3f}")
    
    def validate_nonredundancy(
        self,
        celltype_key: str,
        condition_key: Optional[str] = None,
        count_key: Optional[str] = None,
        n_perm: int = 1000,
        n_matched: int = 100,
        n_perm_matched: int = 200,
        n_boot: int = 1000,
        n_bins: int = 5,
        seed: int = 0,
        verbose: bool = True,
        concordance_key: Optional[str] = "regime_concordance",
        concordance_threshold: float = 0.5,
    ) -> Dict:
        """
        Validate that operator regimes are not redundant with cell types.

        Runs two layers:

        1. Descriptive (backwards-compatible): per-cell-type regime
           diversity and regime × condition fraction tables. Same shape
           existing callers see today.

        2. Statistical (Task 3.1): chi-square + Cramér's V per cell type,
           permutation p-value, Benjamini-Hochberg FDR across cell types,
           bootstrap CIs on regime fractions, AND a depth-matched
           permutation distribution — the acceptance-critical step that
           separates a real condition effect from a UMI-depth artifact.

        Args:
            celltype_key: Key in adata.obs for cell type labels.
            condition_key: Optional key in adata.obs for condition labels.
                Required for the statistical layer.
            count_key: Optional obs column for per-cell UMI count. If None,
                counts are auto-detected via infer_counts_from_adata
                (adata.obs['n_counts'] → 'total_counts' → raw.X sum → X sum).
                Depth-matching is skipped with a warning only if nothing
                is available.
            n_perm: Permutations per cell type for the full-data test.
            n_matched: Random depth-matched subsamples for the matched
                permutation distribution.
            n_perm_matched: Permutations per matched draw.
            n_boot: Bootstrap iterations for regime-fraction CIs.
            n_bins: Depth-quantile bins for depth-matched subsampling.
            seed: Master RNG seed for reproducibility.
            verbose: Print per-cell-type summary as we go.

        Returns:
            Dict with keys:
                celltype_key, condition_key, celltype_regime_diversity,
                condition_comparison (existing backwards-compat contents),
                plus 'statistics' populated when condition_key is provided.
                'statistics' has: per_celltype (chi2, cramers_v, p_perm,
                entropy_mm, n_cells), q_values (BH-corrected), matched
                (depth-matched permutation distribution per cell type),
                bootstrap_ci, depth_matching_applied, n_perm/n_matched/n_boot,
                count_source (str describing where counts came from).
        """
        if celltype_key not in self.adata.obs:
            raise ValueError(f"Cell type key '{celltype_key}' not found in adata.obs")

        celltype = self.adata.obs[celltype_key].values

        validation: Dict = {
            'celltype_key': celltype_key,
            'condition_key': condition_key,
        }

        # ---- Layer 1: descriptive (backwards-compatible) ----
        if verbose:
            print("\n=== Non-Redundancy Test 1: Regime Diversity Within Cell Types ===")
        unique_celltypes = np.unique(celltype)
        celltype_regime_diversity = {}
        for ct in unique_celltypes:
            ct_mask = celltype == ct
            ct_regimes = self.regimes[ct_mask]
            uniq, counts = np.unique(ct_regimes, return_counts=True)
            celltype_regime_diversity[ct] = {
                'n_cells': int(ct_mask.sum()),
                'n_regimes': int(len(uniq)),
                'regime_distribution': dict(zip(uniq, counts.astype(int))),
                # Miller-Madow bias-corrected entropy (Task 3.1). Plug-in
                # entropy was downward-biased at small n, which meant
                # cell-type-level entropy drifted with sample count.
                'entropy': _stats.miller_madow_entropy(counts),
            }
            if verbose:
                print(f"{ct:20s}: {len(uniq)} regimes across {int(ct_mask.sum())} cells")
        validation['celltype_regime_diversity'] = celltype_regime_diversity

        if condition_key is None or condition_key not in self.adata.obs:
            return validation

        conditions = self.adata.obs[condition_key].values
        comparison = self.classifier.compare_regimes_across_conditions(
            self.regimes, conditions, celltype
        )
        validation['condition_comparison'] = comparison

        if verbose:
            print("\n=== Non-Redundancy Test 2: Same Cell Type, Different Conditions ===")
            if 'by_celltype' in comparison:
                for ct in unique_celltypes:
                    if ct in comparison['by_celltype']:
                        print(f"\n{ct}:")
                        for cond, cstats in comparison['by_celltype'][ct].items():
                            fracs = cstats['regime_fractions']
                            print(f"  {cond}: stable={fracs.get('stable', 0):.2f}, "
                                  f"plastic={fracs.get('plastic', 0):.2f}, "
                                  f"unstable={fracs.get('unstable', 0):.2f}")

        # ---- Layer 2: statistical (Task 3.1) ----
        if count_key is not None:
            if count_key not in self.adata.obs.columns:
                raise ValueError(
                    f"count_key='{count_key}' not found in adata.obs.columns"
                )
            counts_arr = np.asarray(self.adata.obs[count_key].values, dtype=float)
            count_source = f"obs['{count_key}']"
        else:
            counts_arr, count_source = _stats.infer_counts_from_adata(self.adata)

        if verbose:
            print("\n=== Non-Redundancy Test 3: Statistical inference (Task 3.1) ===")
            print(f"  count source: {count_source}, "
                  f"perm={n_perm} matched={n_matched} boot={n_boot}")

        result = _stats.test_nonredundancy(
            regimes=self.regimes,
            celltype=celltype,
            condition=conditions,
            counts=counts_arr,
            n_perm=n_perm,
            n_matched=n_matched,
            n_perm_matched=n_perm_matched,
            n_boot=n_boot,
            n_bins=n_bins,
            rng=np.random.default_rng(seed),
        )

        stats_out = result.to_dict()
        stats_out['count_source'] = count_source
        validation['statistics'] = stats_out

        if verbose:
            for ct in unique_celltypes:
                pc = result.per_celltype.get(ct, {})
                q = result.q_values.get(ct, float('nan'))
                mtd = result.matched.get(ct, {})
                print(f"  {ct:20s}: Cramér's V = {pc.get('cramers_v', float('nan')):.3f}, "
                      f"p_perm = {pc.get('p_perm', float('nan')):.4f}, "
                      f"q = {q:.4f}, "
                      f"matched p median = {mtd.get('p_perm_median', float('nan')):.4f}")
            if not result.depth_matching_applied:
                print("  WARNING: no per-cell counts found; depth-matching "
                      "was NOT applied. Effects may be UMI-depth artifacts.")

        # ---- Layer 3: concordance-restricted rerun (Task 3.2) ----
        # If an ensemble has written per-cell concordance to adata.obs,
        # rerun the inference on concordant-only cells. Report BOTH
        # results — the discordant slice is a real cohort and hiding it
        # would misrepresent the ensemble's own uncertainty. See
        # OperatorEnsemble.
        if concordance_key and concordance_key in self.adata.obs.columns:
            conc = np.asarray(self.adata.obs[concordance_key].values, dtype=float)
            keep = conc >= concordance_threshold
            n_keep, n_total = int(keep.sum()), int(len(keep))
            if verbose:
                print(f"\n=== Non-Redundancy Test 4: concordant-only rerun "
                      f"(concordance ≥ {concordance_threshold}) ===")
                print(f"  keeping {n_keep}/{n_total} cells "
                      f"({100 * n_keep / max(n_total, 1):.1f}%)")
            if n_keep >= 8 and len(np.unique(conditions[keep])) >= 2:
                conc_result = _stats.test_nonredundancy(
                    regimes=self.regimes[keep],
                    celltype=celltype[keep],
                    condition=conditions[keep],
                    counts=None if counts_arr is None else counts_arr[keep],
                    n_perm=n_perm,
                    n_matched=n_matched,
                    n_perm_matched=n_perm_matched,
                    n_boot=n_boot,
                    n_bins=n_bins,
                    rng=np.random.default_rng(seed + 1),
                )
                conc_stats = conc_result.to_dict()
                conc_stats['count_source'] = count_source
                conc_stats['concordance_key'] = concordance_key
                conc_stats['concordance_threshold'] = concordance_threshold
                conc_stats['n_cells_kept'] = n_keep
                conc_stats['n_cells_total'] = n_total
                validation['statistics_concordant_only'] = conc_stats
                if verbose:
                    for ct in np.unique(celltype[keep]):
                        pc = conc_result.per_celltype.get(ct, {})
                        q = conc_result.q_values.get(ct, float('nan'))
                        print(f"  {ct:20s}: Cramér's V = {pc.get('cramers_v', float('nan')):.3f}, "
                              f"p_perm = {pc.get('p_perm', float('nan')):.4f}, "
                              f"q = {q:.4f}")
            elif verbose:
                print("  too few cells / conditions after threshold; skipping.")

        return validation
    
    def compare_conditions(
        self,
        condition_key: str,
        celltype_key: Optional[str] = None
    ) -> Dict:
        """
        Compare operator regimes across experimental conditions.
        
        Args:
            condition_key: Key in adata.obs for condition labels
            celltype_key: Optional key for cell type labels
            
        Returns:
            Dictionary with comparison results
        """
        if condition_key not in self.adata.obs:
            raise ValueError(f"Condition key '{condition_key}' not found in adata.obs")
        
        conditions = self.adata.obs[condition_key].values
        celltype = self.adata.obs[celltype_key].values if celltype_key else None
        
        comparison = self.classifier.compare_regimes_across_conditions(
            self.regimes, conditions, celltype
        )
        
        return comparison
    
    def get_regime_statistics(self) -> Dict:
        """Get summary statistics for each operator regime."""
        return self.classifier.get_regime_statistics(self.regimes, self.metrics)
    
    def to_anndata(self) -> ad.AnnData:
        """
        Return AnnData object with atlas results.
        
        Returns:
            AnnData object with operator metrics and regimes stored
        """
        return self.adata
    
    def save(self, filename: str):
        """
        Save atlas results to H5AD file.
        
        Args:
            filename: Output filename (should end with .h5ad)
        """
        self.adata.write_h5ad(filename)
        print(f"Atlas saved to {filename}")
