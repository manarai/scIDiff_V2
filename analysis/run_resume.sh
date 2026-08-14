#!/usr/bin/env bash
# analysis/run_resume.sh — resume run_all.sh from fig2c after path-symlink fix.
#
# Skips items 1-9 (already completed successfully in the initial run):
#   fa_composite_reanalysis.py, init_at_gt_probe.py, top_hvg_null.py,
#   tv_robustness_test.py, concurrent_fix_tests.py, saddle_readouts.py,
#   saddle_dense_sampling.py, fig2b_peak_timing_null_analysis.ipynb,
#   fig2a_synthetic_benchmark.ipynb
#
# Runs items 10-27 (fig2c onward). Appends to the same run_all.status file so
# the full run history is contiguous.

set -u
cd "$(dirname "$0")"

JUPYTER=${JUPYTER:-~/miniconda3/bin/jupyter}
PYTHON=${PYTHON:-~/miniconda3/bin/python}
LOGDIR=logs
STATUS_FILE=run_all.status
mkdir -p "$LOGDIR"

# Remaining items in original run_all.sh order (light → heavy).
declare -a ITEMS=(
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
  "py:supp_bandwidth_sweep/bandwidth_sweep_dc_met.py"
  "nb:figure2_synthetic_and_benchmarks/fig2f_fa_pca_scvi_benchmark.ipynb"
  "nb:figure5_reprogramming/figure5_reprogramming.ipynb"
  "nb:figure5_reprogramming/figure5_reprogramming_tv01.ipynb"
  "py:figure5_reprogramming/make_supp_fig_tv_robustness.py"
  "py:supp_bandwidth_sweep/make_supp_fig_bandwidth.py"
)

echo "=== analysis/run_resume.sh  from fig2c  $(date) ===" | tee -a "$STATUS_FILE"

for entry in "${ITEMS[@]}"; do
  kind="${entry%%:*}"
  path="${entry#*:}"
  short=$(basename "$path")
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

echo "=== RESUME DONE $(date) ===" | tee -a "$STATUS_FILE"
