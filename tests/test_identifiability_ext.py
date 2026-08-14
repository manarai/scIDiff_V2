"""
Tests for the identifiability.py extensions:
  - instability_peak_overlap ground_truth arg;
  - prior_sensitivity_report accuracy vs precision reporting.
"""
from __future__ import annotations

import numpy as np
import torch

from scjdo.validation.identifiability import (
    instability_peak_overlap,
    prior_sensitivity_report,
)


class TestInstabilityPeakGroundTruth:
    def test_consistent_wrong_peak_scores_low_accuracy(self):
        # All runs peak at t=0.9, but ground truth is t=0.5.
        curves = []
        for _ in range(4):
            c = np.zeros(100)
            c[90] = 1.0
            curves.append(c)
        res = instability_peak_overlap(
            curves, peak_window=0.1, ground_truth_peak=0.5
        )
        assert res["overlap_frac"] == 1.0, "Runs agree (precision high)"
        assert res["gt_accuracy"] == 0.0, "But all wrong vs ground truth"
        assert res["gt_mean_error"] > 0.35

    def test_correct_peak_scores_high_accuracy(self):
        curves = []
        for _ in range(3):
            c = np.zeros(100)
            c[50] = 1.0
            curves.append(c)
        res = instability_peak_overlap(
            curves, peak_window=0.05, ground_truth_peak=0.5
        )
        assert res["gt_accuracy"] == 1.0
        assert res["gt_mean_error"] < 0.02

    def test_ground_truth_optional_backcompat(self):
        curves = [np.random.rand(30) for _ in range(3)]
        res = instability_peak_overlap(curves)
        assert "gt_accuracy" not in res
        assert "peak_locations" in res


class TestPriorSensitivityReport:
    def _make_results(self, priors, curves):
        """Helper: attach curves and priors to configs."""
        assert len(priors) == len(curves)
        results = {}
        for i, (prior, curve) in enumerate(zip(priors, curves)):
            results[f"cfg{i}"] = {
                "prior": prior,
                "archetypes": torch.randn(3, 4, 4),
                "instability_curve": curve,
                "regime_fractions": {"stable": 0.5, "plastic": 0.3, "unstable": 0.2},
            }
        return results

    def test_flags_unvaried_prior(self):
        priors = [{"sigma": 0.1}, {"sigma": 0.1}]
        curves = [np.random.rand(30), np.random.rand(30)]
        report = prior_sensitivity_report(self._make_results(priors, curves))
        assert "no listed prior knob varies" in report.lower()

    def test_reports_varied_knob(self):
        priors = [{"sigma": 0.1}, {"sigma": 0.3}, {"sigma": 0.5}]
        c = np.zeros(50)
        c[25] = 1.0
        curves = [c.copy(), c.copy(), c.copy()]
        report = prior_sensitivity_report(self._make_results(priors, curves))
        assert "sigma" in report.lower()
        assert "IDENTIFIABLE" in report

    def test_ground_truth_flows_through(self):
        priors = [{"sigma": s} for s in (0.1, 0.3, 0.5)]
        # All peak at wrong location.
        c = np.zeros(100)
        c[90] = 1.0
        curves = [c.copy() for _ in range(3)]
        report = prior_sensitivity_report(
            self._make_results(priors, curves),
            ground_truth_peak=0.5,
        )
        assert "Ground-truth accuracy" in report
        assert "NOT IDENTIFIABLE" in report

    def test_needs_two_configs(self):
        priors = [{"sigma": 0.1}]
        curves = [np.random.rand(30)]
        report = prior_sensitivity_report(self._make_results(priors, curves))
        assert "need >= 2" in report
