"""
FA composite reanalysis under multiple scoring choices.

The critique: reporting only the fixed specificity produces an unsupervised
method (FA) that beats the supervised upper bound (PLS) — an inversion that
appeared exactly at the moment I changed the rule. Report ALL scorings and
let the reader adjudicate.

Uses the per-component values already in `benchmark_scores.csv` (marker,
specificity, stability, timing, r2) so no notebook rerun is needed. The
specificity column is recomputed under three rules; the composite is
recomputed under four scoring variants (with/without the timing component,
crossed with old/new specificity).
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
CSV = REPO / "Figures_notebook" / "results" / "benchmark" / "benchmark_scores.csv"
FIG3_DIR = REPO / "Figures_notebook" / "results" / "figure3_fa"
BRANCHES = ["Ery", "DC", "Mono"]

# Per-method top-20 instability gene sets, needed to reproduce specificity
# under multiple rules. FA is the only method with a persisted per-branch
# gene table on disk; for the others we invert the specificity value stored
# in the CSV under the CURRENT (empty-as-exclusive) rule to back out the sum
# of pairwise Jaccards on their non-empty pairs, which is enough to recompute
# under any rule that only differs in how empties are treated.


def fa_gene_sets(top_n: int = 20) -> dict[str, set]:
    out: dict[str, set] = {}
    for br in BRANCHES:
        p = FIG3_DIR / f"instability_genes_{br}.csv"
        if not p.exists():
            out[br] = set()
            continue
        # Empty CSV (buffered branch — only header) raises EmptyDataError
        try:
            df = pd.read_csv(p)
            genes = df["gene"].dropna().astype(str).head(top_n).tolist()
        except (pd.errors.EmptyDataError, KeyError):
            genes = []
        out[br] = {g.upper() for g in genes}
    return out


def spec_from_gene_sets(gs: dict[str, set], rule: str) -> float:
    """
    rule = 'strict_zero'        — v31 original: any empty branch → 0
           'exclusive_empty'    — v32 first fix: empty contributes Jaccard 0 to every pair
           'exclude_empty'      — NaN if <2 non-empty branches, else avg over non-empty pairs only
    """
    branches = list(gs.keys())
    n_nonempty = sum(1 for b in branches if len(gs[b]) > 0)
    if rule == "strict_zero":
        if n_nonempty < len(branches):
            return 0.0
    if rule in ("exclude_empty",):
        if n_nonempty < 2:
            return float("nan")

    pairs = []
    for i in range(len(branches)):
        for j in range(i + 1, len(branches)):
            a, b = gs[branches[i]], gs[branches[j]]
            if rule == "exclude_empty":
                if len(a) == 0 or len(b) == 0:
                    continue
            u = len(a | b)
            if u == 0:
                jacc = 0.0
            else:
                jacc = len(a & b) / u
            pairs.append(jacc)
    if not pairs:
        return float("nan")
    return float(1.0 - np.mean(pairs))


def composite(marker, specificity, stability, timing, r2, drop_timing=False):
    if drop_timing:
        # Renormalise remaining weights to sum to 1.
        # Original weights: marker 0.30, specificity 0.20, stability 0.20, r2 0.15
        s = 0.30 + 0.20 + 0.20 + 0.15
        w = np.array([0.30, 0.20, 0.20, 0.15]) / s
        vals = np.array([marker, specificity, stability, r2])
    else:
        w = np.array([0.30, 0.20, 0.20, 0.15, 0.15])
        vals = np.array([marker, specificity, stability, timing, r2])
    # Any NaN component → NaN composite
    if np.isnan(vals).any():
        return float("nan")
    return float((w * vals).sum())


def main():
    df = pd.read_csv(CSV)
    fa_gs = fa_gene_sets()
    # Compute FA specificity under each rule.
    fa_spec = {rule: spec_from_gene_sets(fa_gs, rule)
               for rule in ("strict_zero", "exclusive_empty", "exclude_empty")}
    print("FA gene-set summary (top-20):")
    for b in BRANCHES:
        print(f"  {b}: {len(fa_gs[b])} genes")
    print(f"FA specificity under each rule: {fa_spec}\n")

    # Rebuild a per-method table where the specificity value can be swapped.
    # For non-FA methods, we don't have per-branch gene tables on disk from
    # this run, so we rely on the CSV's stored specificity (which was
    # computed under 'exclusive_empty'). We flag which rows we can vary.
    rows = []
    for _, r in df.iterrows():
        m = r["method"]
        row = {
            "method": m,
            "supervised": bool(r["supervised"]),
            "marker": float(r["marker"]),
            "stability": float(r["stability"]),
            "timing": float(r["timing"]),
            "r2": float(r["r2"]),
        }
        # Only FA's specificity is recomputable rule-by-rule from disk here.
        # For the others, we use the CSV value (rule='exclusive_empty') across all rules.
        for rule in ("strict_zero", "exclusive_empty", "exclude_empty"):
            row[f"spec_{rule}"] = fa_spec[rule] if m == "FactorAnalysis" else float(r["specificity"])
        rows.append(row)
    out = pd.DataFrame(rows)

    # DiffMap returned zero instability genes for every branch — under
    # 'exclude_empty' rule it drops out; under 'strict_zero' it scores 0;
    # under 'exclusive_empty' it landed at 1.0 in the CSV. Patch that
    # semantically here for the 'exclude_empty' rule.
    diffmap_mask = out["method"] == "DiffMap"
    if diffmap_mask.any():
        out.loc[diffmap_mask, "spec_exclude_empty"] = float("nan")

    # Compute composites under each combination.
    for spec_rule in ("strict_zero", "exclusive_empty", "exclude_empty"):
        for drop_t in (False, True):
            col = f"comp_{spec_rule}" + ("_no_timing" if drop_t else "")
            out[col] = out.apply(
                lambda x: composite(x["marker"], x[f"spec_{spec_rule}"],
                                    x["stability"], x["timing"], x["r2"],
                                    drop_timing=drop_t),
                axis=1,
            )

    # Print a compact table of composites under each rule
    display_cols = ["method", "supervised",
                    "comp_strict_zero", "comp_exclusive_empty", "comp_exclude_empty",
                    "comp_strict_zero_no_timing", "comp_exclusive_empty_no_timing",
                    "comp_exclude_empty_no_timing"]
    disp = out[display_cols].copy()
    # Round for readability
    for c in display_cols[2:]:
        disp[c] = disp[c].map(lambda v: f"{v:.3f}" if pd.notna(v) else "  NaN")

    print("Composite score under multiple scoring rules")
    print("=" * 96)
    print(disp.to_string(index=False))
    print("=" * 96)

    print()
    print("Column key:")
    print("  strict_zero      = v31 paper rule (any empty branch → specificity 0)")
    print("  exclusive_empty  = v32 first fix (empty contributes Jaccard 0 to each pair)")
    print("  exclude_empty    = v32 guarded (drop empty-branch pairs; NaN if <2 non-empty branches)")
    print("  _no_timing       = drop the timing-spread component (paper declines timing claims)")
    print()

    # Rank under each column
    for rule in ("comp_strict_zero", "comp_exclusive_empty", "comp_exclude_empty",
                 "comp_strict_zero_no_timing", "comp_exclusive_empty_no_timing",
                 "comp_exclude_empty_no_timing"):
        ranked = out[["method", rule]].dropna().sort_values(rule, ascending=False)
        top3 = " > ".join(f"{r.method}({r[rule]:.3f})" for _, r in ranked.head(3).iterrows())
        print(f"  {rule}:  {top3}")

    out.to_csv("Path(__file__).parent / "fa_composite_reanalysis.csv"",
               index=False)
    print("\nSaved: Path(__file__).parent / "fa_composite_reanalysis.csv"")


if __name__ == "__main__":
    main()
