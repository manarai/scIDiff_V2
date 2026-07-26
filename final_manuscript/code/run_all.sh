#!/usr/bin/env bash
# analysis/run_all.sh — sequential rerunner for all manuscript-figure
# notebooks and diagnostic scripts.
#
# Ordering: light → heavy. Total wall-clock is ~15-20 hours; the two heaviest
# notebooks are fig5_reprogramming and fig5_reprogramming_tv01 (~7 hrs each)
# and fig2f_fa_pca_scvi_benchmark (~8.5 hrs). Per-item status lands in
# analysis/run_all.status; per-item logs in analysis/logs/.
#
# Usage:
#   analysis/run_all.sh             # run everything, sequentially
#   analysis/run_all.sh fast        # skip the three heaviest notebooks
#   analysis/run_all.sh figN         # run only the items for one figure

set -u
cd "$(dirname "$0")"

JUPYTER=${JUPYTER:-~/miniconda3/bin/jupyter}
PYTHON=${PYTHON:-~/miniconda3/bin/python}
LOGDIR=logs
STATUS_FILE=run_all.status
mkdir -p "$LOGDIR"
: > "$STATUS_FILE"

MODE="${1:-all}"

# Sequential list: (kind, path). kind = "nb" or "py".
# Order: fastest first so failures on cheap items surface quickly.
declare -a ITEMS=(
  # Fast diagnostics (< 1 min each)
  "py:figure2_synthetic_and_benchmarks/fa_composite_reanalysis.py"
  "py:figure2_synthetic_and_benchmarks/init_at_gt_probe.py"
  "py:figure4_concordance/top_hvg_null.py"

  # Fast Python experiments (a few min each)
  "py:figure2_synthetic_and_benchmarks/tv_robustness_test.py"
  "py:figure2_synthetic_and_benchmarks/concurrent_fix_tests.py"
  "py:supp_saddle_localization/saddle_readouts.py"
  "py:supp_saddle_localization/saddle_dense_sampling.py"

  # Fast notebook (2-3 min)
  "nb:figure2_synthetic_and_benchmarks/fig2b_peak_timing_null_analysis.ipynb"

  # Medium notebooks (10-40 min each)
  "nb:figure2_synthetic_and_benchmarks/fig2a_synthetic_benchmark.ipynb"
  "nb:figure2_synthetic_and_benchmarks/fig2c_embedding_benchmark.ipynb"
  "nb:figure2_synthetic_and_benchmarks/fig2d_pseudotime_embedding_benchmark.ipynb"
  "nb:figure2_synthetic_and_benchmarks/fig2e_fa_benchmark.ipynb"
  "nb:figure3_hematopoiesis/figure3_hematopoiesis.ipynb"
  "nb:figure4_concordance/fig4_01_scjdo_operators.ipynb"
  "nb:figure4_concordance/fig4_02_splicejac_operators.ipynb"
  "nb:figure4_concordance/fig4_03_dynamo_operators.ipynb"
  "nb:figure4_concordance/fig4_04_concordance_metrics.ipynb"
  "nb:figure6_k562_perturbseq/figure6_k562.ipynb"
  "nb:supp_saddle_localization/saddle_localization_v3_label_blind.ipynb"
  "nb:supp_multiome_note/supp_multiome_FA_integration.ipynb"
  "nb:supp_multiome_note/supp_multiome_drift_cellcycle_regressed.ipynb"

  # Bandwidth-sweep (also medium — 15 min)
  "py:supp_bandwidth_sweep/bandwidth_sweep_dc_met.py"

  # HEAVY notebooks (hours each)
  "nb:figure2_synthetic_and_benchmarks/fig2f_fa_pca_scvi_benchmark.ipynb"
  "nb:figure5_reprogramming/figure5_reprogramming.ipynb"
  "nb:figure5_reprogramming/figure5_reprogramming_tv01.ipynb"

  # Post-processing figure builders (depend on outputs above)
  "py:figure5_reprogramming/make_supp_fig_tv_robustness.py"
  "py:supp_bandwidth_sweep/make_supp_fig_bandwidth.py"
)

HEAVY_ITEMS=(
  "nb:figure2_synthetic_and_benchmarks/fig2f_fa_pca_scvi_benchmark.ipynb"
  "nb:figure5_reprogramming/figure5_reprogramming.ipynb"
  "nb:figure5_reprogramming/figure5_reprogramming_tv01.ipynb"
)

is_heavy() {
  local item="$1"
  for h in "${HEAVY_ITEMS[@]}"; do
    [[ "$h" == "$item" ]] && return 0
  done
  return 1
}

should_run() {
  local item="$1"
  case "$MODE" in
    all)   return 0 ;;
    fast)  is_heavy "$item" && return 1 || return 0 ;;
    fig*)  [[ "$item" == *"$MODE"* ]] && return 0 || return 1 ;;
    *)     return 0 ;;
  esac
}

echo "=== analysis/run_all.sh  MODE=$MODE  $(date) ===" | tee -a "$STATUS_FILE"

for entry in "${ITEMS[@]}"; do
  kind="${entry%%:*}"
  path="${entry#*:}"
  short=$(basename "$path")
  if ! should_run "$entry"; then
    echo "SKIP  $short  (mode=$MODE)" | tee -a "$STATUS_FILE"
    continue
  fi
  log="$LOGDIR/${short%.*}.log"
  t0=$SECONDS
  echo "START $short  $(date +%H:%M:%S)" | tee -a "$STATUS_FILE"
  if [[ "$kind" == "nb" ]]; then
    "$JUPYTER" nbconvert --to notebook --execute "$path" \
      --output "$(basename "$path")" \
      --ExecutePreprocessor.timeout=32400 \
      >"$log" 2>&1
    rc=$?
  else
    # Python script — run from its own directory so relative paths work
    d=$(dirname "$path")
    f=$(basename "$path")
    (cd "$d" && "$PYTHON" -u "$f") >"$log" 2>&1
    rc=$?
  fi
  dt=$((SECONDS - t0))
  if [[ $rc -eq 0 ]]; then
    echo "OK    $short  ${dt}s" | tee -a "$STATUS_FILE"
  else
    echo "FAIL  $short  ${dt}s  rc=$rc  (see $log)" | tee -a "$STATUS_FILE"
  fi
done

echo "=== DONE $(date) ===" | tee -a "$STATUS_FILE"
