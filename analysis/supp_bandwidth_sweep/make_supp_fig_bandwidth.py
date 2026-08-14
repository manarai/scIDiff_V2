"""
Supplementary Fig. Sxxx-BW — bandwidth-sweep sensitivity for the MET peak.

Renders the SCP295 MET-branch bandwidth sweep as a two-panel figure:
- Left: peak τ vs h, with the manuscript's reported peak τ overlaid
- Right: peak Re(λ_max) vs h (crosses zero — most bandwidths return negative peak)
Both panels annotate n_eff at each supported peak.
"""
from __future__ import annotations
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

REPO = Path("/Users/terooatt/Downloads/scJDO")
d = json.loads((REPO / "scratchpad_run/bandwidth_sweep_dc_met.json").read_text())
met = d["MET"]

mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.titlesize": 10, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
})

hs = [r["bandwidth"] for r in met]
peak_ts = [r["peak_t"] for r in met]
peak_lams = [r["peak_lam"] for r in met]
n_effs = [r["n_eff_at_peak"] for r in met]

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

# Panel A — peak τ
ax = axes[0]
ax.plot(hs, peak_ts, "o-", color="#355C9A", ms=7, lw=1.5, label="peak τ (this sweep)")
ax.axhline(0.768, color="#E8884A", ls="--", alpha=0.7,
           label="manuscript Fig 5 peak τ = 0.77")
for h, t, n in zip(hs, peak_ts, n_effs):
    ax.annotate(f"n_eff={n:.0f}", (h, t), textcoords="offset points",
                xytext=(6, 6), fontsize=7, color="#666666")
ax.set_xlabel("Kernel bandwidth h")
ax.set_ylabel("Peak τ (interior, n_eff ≥ 20 filter)")
ax.set_xscale("log")
ax.set_ylim(-0.05, 1.05)
ax.set_title("A | Peak τ vs bandwidth", loc="left", fontweight="bold")
ax.grid(True, alpha=0.3)
ax.legend(loc="lower right", fontsize=8)

# Panel B — peak Re(λ_max)
ax = axes[1]
ax.plot(hs, peak_lams, "s-", color="#2A9D8F", ms=7, lw=1.5, label="peak Re(λ_max)")
ax.axhline(0, color="black", ls=":", alpha=0.5, label="λ = 0")
ax.axhline(0.05, color="red", ls=":", alpha=0.5,
           label="v31 'locally sensitive' threshold")
for h, l, n in zip(hs, peak_lams, n_effs):
    ax.annotate(f"n_eff={n:.0f}", (h, l), textcoords="offset points",
                xytext=(6, 6), fontsize=7, color="#666666")
ax.set_xlabel("Kernel bandwidth h")
ax.set_ylabel("Peak Re(λ_max)")
ax.set_xscale("log")
ax.set_title("B | Peak Re(λ_max) vs bandwidth\n(only h ≈ 0.01 gives a positive peak eigenvalue)",
             loc="left", fontweight="bold")
ax.grid(True, alpha=0.3)
ax.legend(loc="upper right", fontsize=8)

fig.suptitle(
    "Supplementary Fig. Sxxx-BW  |  SCP295 MET-branch bandwidth sweep\n"
    "Peak location shifts monotonically 0.60 → 0.95 across h ∈ [0.01, 0.15]; "
    "only h = 0.01 gives a positive peak eigenvalue on this simplified pipeline "
    "(3344-cell MET subset, 1000 epochs). Any timing claim should be checked "
    "against a bandwidth sweep on the exact same pipeline.",
    fontsize=9, y=1.02,
)
plt.tight_layout()
out_pdf = REPO / "scratchpad_run/supp_fig_bandwidth_sweep.pdf"
out_png = REPO / "scratchpad_run/supp_fig_bandwidth_sweep.png"
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, bbox_inches="tight", dpi=200)
print(f"Saved:\n  {out_pdf}\n  {out_png}")
