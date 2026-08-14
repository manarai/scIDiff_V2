"""
Ensemble of drift-field members for gauge-invariant operator regimes.

Task 3.2. Any single trained drift field is one arbitrary representative
from a gauge-admissible set — snapshot data determines the drift only up
to a divergence-free (rotational) circulation, and scJDO resolves that
ambiguity by minimizing kinetic energy. Which is a regularization
choice, not a measurement.

Ensembling N members trained with different seeds / initializations /
prior settings samples that admissible set. Cells whose regime agrees
across the ensemble are approximately gauge-invariant; cells whose
regime flips between members are gauge artifacts. Concordance across an
ensemble is therefore the closest thing to an honest per-cell
uncertainty quantifier the library can offer.

Key outputs written to adata.obs by :meth:`OperatorEnsemble.build`:

- ``ensemble_regime`` — modal regime per cell (majority vote).
- ``regime_concordance`` — fraction of members that assigned the modal
  regime; a real uncertainty quantifier. Replaces the interim
  ``boundary_distance`` (Task 2.3) for callers who have an ensemble.

Global outputs:

- ``discordant_fraction`` — fraction of cells with concordance below a
  configurable threshold. Reported prominently by :meth:`print_summary`.
  If more than roughly a third of cells are discordant, the atlas is not
  measuring a stable property and that finding must be stated — do not
  bury it.

Design note. This class does NOT train members. Training a scJDO drift
field is heavy and typically done outside a Python session; forcing
training into this class would couple the ensemble to a specific
trainer configuration. Accept a list of pre-trained models via
:meth:`__init__`; document what SHOULD vary across members
(seed AND initialization AND velocity-prior setting — seed-only
variation underestimates the gauge-freedom spread).
"""

from __future__ import annotations

from collections import Counter
from typing import Dict, List, Optional, Sequence

import numpy as np
import torch

from .operator_metrics import OperatorMetrics
from .regime_classifier import OperatorRegimeClassifier


REGIME_NAMES = ("stable", "plastic", "unstable", "deeply_stable")


class OperatorEnsemble:
    """
    Aggregate operator regimes across N drift-model members.

    Args:
        drift_models: List of pre-trained drift models. Each must expose
            ``forward(x, t)`` returning the drift at ``(x, t)`` — the
            same contract :class:`OperatorMetrics` uses.
        epsilon: Threshold for near-neutral modes (passed to each
            member's :class:`OperatorMetrics`).
        device: 'cpu' / 'cuda' / etc.
        symmetric_only: If True, each member eigendecomposes the
            symmetric part (J + J.T)/2. See :class:`OperatorMetrics`.
        stability_quantile: Per-cell quantile of Re(λ) used as the
            dimension-stable stability metric. See :class:`OperatorMetrics`.

    What to vary across members (a note for the person training them):
        - Random seed (necessary but not sufficient).
        - Network initialization (implicit in seed unless you fix seeds
          and vary init explicitly).
        - Velocity-prior ON vs OFF (single most important axis; controls
          how the model resolves the gauge freedom).
        - Bandwidth / kernel settings if you use adaptive windowing.
        Members that differ only by seed underestimate the true spread.
    """

    def __init__(
        self,
        drift_models: Sequence,
        epsilon: float = 0.1,
        device: str = "cpu",
        symmetric_only: bool = False,
        stability_quantile: float = 0.05,
    ):
        if len(drift_models) < 2:
            raise ValueError(
                "OperatorEnsemble requires at least 2 members; got "
                f"{len(drift_models)}. A one-member 'ensemble' cannot "
                "quantify gauge ambiguity."
            )
        self.drift_models = list(drift_models)
        self.epsilon = epsilon
        self.device = device
        self.symmetric_only = symmetric_only
        self.stability_quantile = stability_quantile

        # Populated by build().
        self.per_member_regimes: Optional[np.ndarray] = None  # (n_members, n_cells)
        self.per_member_metrics: Optional[List[Dict[str, np.ndarray]]] = None
        self.ensemble_regime: Optional[np.ndarray] = None     # (n_cells,)
        self.concordance: Optional[np.ndarray] = None         # (n_cells,)
        self.adata = None

    @property
    def n_members(self) -> int:
        return len(self.drift_models)

    def build(
        self,
        adata,
        use_rep: str = "X_pca",
        pseudotime_key: str = "pseudotime",
        batch_size: int = 32,
        thresholds: str = "hyperbolic",
        threshold_unstable: Optional[float] = None,
        threshold_plastic: Optional[float] = None,
        threshold_deeply_stable: Optional[float] = None,
        plasticity_threshold: float = 0.3,
        deep_stable_scale_k: float = 2.0,
        store_member_metrics: bool = False,
        write_to_obs: bool = True,
        verbose: bool = True,
    ) -> None:
        """
        Compute regimes for every member and aggregate.

        Args:
            adata: AnnData with ``obsm[use_rep]`` and ``obs[pseudotime_key]``.
            use_rep, pseudotime_key: see :class:`StableOperatorAtlas`.
            batch_size, thresholds, threshold_*: passed through to each
                member's metric computer / classifier. Every member uses
                the SAME thresholds so concordance is a pure ensemble-
                gauge measurement, not a threshold-noise measurement.
            store_member_metrics: If True, keep per-member metric dicts
                on ``self.per_member_metrics``. Memory-heavy for large N;
                off by default.
            write_to_obs: If True, mirror the ensemble regime and
                concordance to ``adata.obs``.
        """
        if use_rep not in adata.obsm:
            raise ValueError(f"'{use_rep}' not in adata.obsm")
        if pseudotime_key not in adata.obs.columns:
            raise ValueError(
                f"pseudotime_key '{pseudotime_key}' not in adata.obs — "
                "scOpAtlas requires a real pseudotime axis (see "
                "StableOperatorAtlas._validate_inputs)."
            )

        self.adata = adata
        X = torch.tensor(adata.obsm[use_rep], dtype=torch.float32)
        pseudotime = adata.obs[pseudotime_key].values

        per_member_regimes = []
        per_member_metrics = [] if store_member_metrics else None

        for k, model in enumerate(self.drift_models):
            if verbose:
                print(f"[ensemble] member {k + 1}/{self.n_members} …", flush=True)
            metrics_computer = OperatorMetrics(
                drift_model=model,
                epsilon=self.epsilon,
                device=self.device,
                symmetric_only=self.symmetric_only,
                stability_quantile=self.stability_quantile,
            )
            m = metrics_computer.compute_metrics_at_pseudotime(
                X, pseudotime, batch_size=batch_size,
            )
            classifier = OperatorRegimeClassifier(
                thresholds=thresholds,
                threshold_unstable=threshold_unstable,
                threshold_plastic=threshold_plastic,
                threshold_deeply_stable=threshold_deeply_stable,
                plasticity_threshold=plasticity_threshold,
                deep_stable_scale_k=deep_stable_scale_k,
            )
            r, _masks = classifier.classify(m)
            per_member_regimes.append(np.asarray(r, dtype=object))
            if store_member_metrics:
                per_member_metrics.append(m)

        self.per_member_regimes = np.stack(per_member_regimes, axis=0)  # (M, N)
        self.per_member_metrics = per_member_metrics

        self.ensemble_regime, self.concordance = _aggregate_modal_and_concordance(
            self.per_member_regimes
        )

        if write_to_obs:
            adata.obs["ensemble_regime"] = self.ensemble_regime
            adata.obs["regime_concordance"] = self.concordance

        if verbose:
            self.print_summary()

    # ------------------------------------------------------------------ #
    # Summaries                                                           #
    # ------------------------------------------------------------------ #

    def discordant_fraction(self, threshold: float = 0.5) -> float:
        """Fraction of cells whose concordance falls below `threshold`.

        Report this number prominently. If it exceeds ~1/3, the atlas is
        not measuring a stable property.
        """
        if self.concordance is None:
            raise RuntimeError("call .build() first")
        return float((self.concordance < threshold).mean())

    def print_summary(self, thresholds: Sequence[float] = (0.5, 0.75, 1.0)) -> None:
        """Print modal-regime counts + discordant-fraction at several
        thresholds. The 1.0 row is the strictest reading ("every member
        agreed"); the 0.5 row is the usual majority-vote reading.
        """
        if self.ensemble_regime is None:
            raise RuntimeError("call .build() first")
        print("\n=== Ensemble Summary ===")
        n = len(self.ensemble_regime)
        print(f"Members: {self.n_members}, cells: {n}")
        print("\nModal regime counts (majority vote across members):")
        uniq, counts = np.unique(self.ensemble_regime, return_counts=True)
        for name, k in zip(uniq, counts):
            print(f"  {name:15s}: {int(k):6d} cells ({k / n * 100:5.1f}%)")
        print("\nDiscordant-cell fractions (unstable across ensemble):")
        for t in thresholds:
            # 'discordant' == concordance strictly below t. At t=1.0 that's
            # 'not fully unanimous'.
            frac = float((self.concordance < t).mean())
            flag = "  <-- flag" if (t == 0.5 and frac > 1 / 3) else ""
            print(f"  concordance < {t:.2f}: {frac * 100:5.1f}%{flag}")

    def concordant_mask(self, threshold: float = 0.5) -> np.ndarray:
        """Boolean mask of cells with concordance ≥ threshold."""
        if self.concordance is None:
            raise RuntimeError("call .build() first")
        return self.concordance >= threshold


def _aggregate_modal_and_concordance(
    per_member_regimes: np.ndarray,
) -> tuple:
    """
    Compute per-cell modal regime and concordance from
    (n_members, n_cells) labels.

    Concordance = count(modal) / n_members, in [1/n_members, 1].

    Ties are broken deterministically by the ordering returned by
    Counter.most_common — which itself falls back to insertion order.
    We insert regimes in REGIME_NAMES order so ties resolve
    predictably across runs. (Ties on 2 out of 3 members with 4 regime
    labels are common in a small ensemble and this deterministic tie-
    break keeps runs reproducible.)
    """
    n_members, n_cells = per_member_regimes.shape
    modal = np.empty(n_cells, dtype=object)
    concordance = np.zeros(n_cells, dtype=float)
    order_map = {name: i for i, name in enumerate(REGIME_NAMES)}
    for i in range(n_cells):
        col = per_member_regimes[:, i]
        c = Counter(col)
        # Sort by (-count, order_map) so highest count wins; ties broken
        # by the canonical regime ordering.
        best = sorted(c.items(), key=lambda kv: (-kv[1], order_map.get(kv[0], 999)))[0]
        modal[i] = best[0]
        concordance[i] = best[1] / n_members
    return modal, concordance
