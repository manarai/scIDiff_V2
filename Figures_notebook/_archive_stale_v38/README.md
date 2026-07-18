# Stale artifacts archived during v38 audit

This directory holds files that were flagged as stale during the manuscript-v38
notebook↔code audit. They are moved here rather than deleted so anything
inadvertently useful can be recovered. Nothing here should be relied on by any
currently-shipped notebook.

## Contents

- `Figure3_FA.ipynb.bak` — May 28 backup of `Figure3_FA.ipynb`. Was the only
  git-tracked stale file (`git mv`'d into place).
- `Temp/` — six obsolete Figure-3 draft notebooks
  (`Figure3.ipynb`, `Figure3_updated.ipynb`, `Figure3_updated_2.ipynb`,
  `Figure3_palantir.ipynb`, `Figure3_FA_benchmark.ipynb`,
  `Figure3_FA-Copy1.ipynb`) plus its own checkpoint dir.
- `run_all.log` — May 31 log from an earlier full-pipeline run.
- `results/figure3/` and `results/figure3_v2/` — outputs keyed to a superseded
  cluster labelling (`W36975`, `W39076`, `W37276`, `W37401`, `W38095`,
  `W38285`). No notebook in the current tree regenerates these.
- `results/figure4_smoke/` — smoke-test outputs superseded by
  `Figures_notebook/results/figure4/`.
- `ipynb_checkpoints/` — three checkpoints whose "live" notebook was renamed
  or removed (`Figure3-checkpoint.ipynb`,
  `Figure3_FA_benchmark-checkpoint.ipynb`,
  `Figure3_FA.ipynb-checkpoint.bak`). Live-notebook checkpoints were left in
  `Figures_notebook/.ipynb_checkpoints/`.
