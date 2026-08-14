"""
scOpAtlas: Stable Operator Atlas for scJDO

This module provides operator-based state definition for single-cell data.
It defines cellular states by their local stability structure rather than
expression patterns alone.

Main components:
- OperatorMetrics: Compute eigenvalue-derived metrics
- OperatorRegimeClassifier: Classify cells into operator regimes
- StableOperatorAtlas: Main atlas object
- OperatorEmbedding: Project operators into low-dimensional space
- OperatorClustering: Clustering utilities for operator space
"""

from .operator_metrics import OperatorMetrics
from .regime_classifier import OperatorRegimeClassifier
from .atlas_builder import StableOperatorAtlas
from .operator_embedding import OperatorEmbedding, compute_operator_embedding
from .clustering import OperatorClustering, quick_operator_clustering
from .flow_topology import compute_ftle, assign_basins, fork_detector
from .ensemble import OperatorEnsemble  # Task 3.2 — ensemble across drift-model members.
from .quality import check_drift_quality, DriftQualityError  # Task 3.3 — drift-quality gate.
from . import statistics  # Task 3.1 inference module — exposed as submodule.
from . import quality  # Task 3.3 quality-gate module — exposed as submodule.

__all__ = [
    "OperatorMetrics",
    "OperatorRegimeClassifier",
    "StableOperatorAtlas",
    "OperatorEmbedding",
    "compute_operator_embedding",
    "OperatorClustering",
    "quick_operator_clustering",
    "compute_ftle",
    "assign_basins",
    "fork_detector",
    "OperatorEnsemble",
    "check_drift_quality",
    "DriftQualityError",
    "statistics",
    "quality",
]
