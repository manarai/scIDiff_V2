"""
scjdo.validation
==================
Reviewer-response validation utilities for the scJDO manuscript.

Submodules
----------
null_models
    Temporal shuffle and continuous-control null models to verify that
    archetype coordination motifs are not artifacts of low-rank
    factorization.

robustness
    Gene-level Jaccard overlap and Spearman rank correlation across
    different pseudotime methods (DPT, Palantir, Slingshot).

identifiability
    Utilities for quantifying what is and is not identifiable in the
    scJDO operator framework: archetype cosine similarity across runs,
    instability peak timing, and model-sensitivity reports.
"""

from scjdo.validation.identifiability import (
    PRIOR_KNOBS,
    archetype_cosine_similarity,
    instability_peak_overlap,
    model_sensitivity_report,
    prior_sensitivity_report,
)
from scjdo.validation.null_models import (
    continuous_control_null,
    run_null_comparison,
    temporal_shuffle_null,
)
from scjdo.validation.robustness import (
    gene_overlap_across_pseudotimes,
    pseudotime_sensitivity_report,
)

__all__ = [
    # null models
    "temporal_shuffle_null",
    "continuous_control_null",
    "run_null_comparison",
    # robustness
    "gene_overlap_across_pseudotimes",
    "pseudotime_sensitivity_report",
    # identifiability
    "archetype_cosine_similarity",
    "instability_peak_overlap",
    "model_sensitivity_report",
    "prior_sensitivity_report",
    "PRIOR_KNOBS",
]
