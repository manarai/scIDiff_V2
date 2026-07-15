"""
Standard figure panels for scJDO drift field analysis.
All functions read from adata.uns[key] populated by sjd.tl.fit_drift.
"""
from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np


_ARCH_COLORS = ["#E63946", "#2A9D8F", "#E9C46A", "#457B9D", "#8338EC",
                "#F4A261", "#264653", "#E76F51", "#A8DADC", "#6D6875"]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _get(adata, key, field):
    if key not in adata.uns:
        raise KeyError(f"'{key}' not in adata.uns. Run sjd.tl.fit_drift first.")
    return adata.uns[key][field]


def _grid_quiver(embed: np.ndarray, drift: np.ndarray, n_bins: int = 20, min_cells: int = 3):
    x, y = embed[:, 0], embed[:, 1]
    xmin, xmax = x.min() - .05, x.max() + .05
    ymin, ymax = y.min() - .05, y.max() + .05
    gx = np.linspace(xmin, xmax, n_bins)
    gy = np.linspace(ymin, ymax, n_bins)
    gu = np.zeros((n_bins, n_bins))
    gv = np.zeros((n_bins, n_bins))
    gn = np.zeros((n_bins, n_bins), dtype=int)
    xi = np.clip(((x - xmin) / (xmax - xmin) * (n_bins - 1)).astype(int), 0, n_bins - 1)
    yi = np.clip(((y - ymin) / (ymax - ymin) * (n_bins - 1)).astype(int), 0, n_bins - 1)
    np.add.at(gu, (yi, xi), drift[:, 0])
    np.add.at(gv, (yi, xi), drift[:, 1])
    np.add.at(gn, (yi, xi), 1)
    mask = gn >= min_cells
    gu[mask] /= gn[mask]
    gv[mask] /= gn[mask]
    GX, GY = np.meshgrid(gx, gy)
    return GX[mask], GY[mask], gu[mask], gv[mask]


# ---------------------------------------------------------------------------
# Panel functions
# ---------------------------------------------------------------------------

def drift_field(
    adata,
    key: str = "scjdo",
    basis: str = "X_pca",
    color: str = "pseudotime",
    velocity_key: Optional[str] = None,
    stream: bool = True,
    stream_density: float = 1.2,
    n_grid: int = 30,
    min_cells: int = 3,
    ax: Optional[plt.Axes] = None,
    save: Optional[str] = None,
    **scatter_kw,
):
    """
    Plot drift field on the specified embedding.

    By default uses **streamlines** (``stream=True``), which integrate the
    flow continuously and make bifurcations visible even when the 2D projection
    captures only a fraction of the total velocity energy.

    Parameters
    ----------
    basis         : 'X_pca' or 'X_umap'.
    color         : 'pseudotime' (default) or any column in adata.obs.
    velocity_key  : Velocity to visualize. Auto-selects 'X_velocity_pseudo'
                    (pseudotime gradient, always coherent) when basis='X_pca',
                    or 'X_drift' (full model drift) otherwise.
    stream        : Use streamplot (True, default) or quiver arrows (False).
    stream_density: Streamline density passed to matplotlib.streamplot.
    n_grid        : Grid resolution for velocity interpolation.
    """
    from scipy.interpolate import griddata

    if basis not in adata.obsm:
        raise KeyError(f"'{basis}' not in adata.obsm.")

    # ── Choose velocity source ─────────────────────────────────────────────
    if velocity_key is None:
        if basis == "X_pca" and "X_velocity_pseudo" in adata.obsm:
            velocity_key = "X_velocity_pseudo"   # always coherent
        elif "X_drift" in adata.obsm:
            velocity_key = "X_drift"
        else:
            raise KeyError("No drift vectors found. Run sjd.tl.fit_drift first.")

    if velocity_key not in adata.obsm:
        raise KeyError(f"'{velocity_key}' not in adata.obsm.")

    embed = adata.obsm[basis][:, :2]
    vel   = adata.obsm[velocity_key][:, :2]

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(5.5, 4.5))

    # ── Scatter cells ──────────────────────────────────────────────────────
    if color == "pseudotime" and "pseudotime" in adata.obs.columns:
        c_vals = adata.obs["pseudotime"].values
        sc = ax.scatter(embed[:, 0], embed[:, 1], c=c_vals, cmap="plasma",
                        s=scatter_kw.pop("s", 8), alpha=scatter_kw.pop("alpha", 0.8),
                        linewidths=0, **scatter_kw)
        plt.colorbar(sc, ax=ax, label="Pseudotime", fraction=0.046, pad=0.04)
    elif color in adata.obs.columns:
        from matplotlib.lines import Line2D
        cats = adata.obs[color].astype(str).values
        uniq = sorted(set(cats))
        cmap = plt.cm.get_cmap("tab20", len(uniq))
        ci   = {c: i for i, c in enumerate(uniq)}
        cols = [cmap(ci[c]) for c in cats]
        ax.scatter(embed[:, 0], embed[:, 1], c=cols,
                   s=scatter_kw.pop("s", 8), alpha=scatter_kw.pop("alpha", 0.75),
                   linewidths=0, **scatter_kw)
        handles = [plt.Line2D([0], [0], marker="o", color="w",
                               markerfacecolor=cmap(ci[c]), markersize=6, label=c)
                   for c in uniq]
        ax.legend(handles=handles, fontsize=6, ncol=2, loc="upper right", framealpha=0.8)
    else:
        ax.scatter(embed[:, 0], embed[:, 1], s=8, alpha=0.7, color="#aaa")

    # ── Drift overlay ──────────────────────────────────────────────────────
    xmin, xmax = embed[:, 0].min(), embed[:, 0].max()
    ymin, ymax = embed[:, 1].min(), embed[:, 1].max()
    pad = 0.05
    gx = np.linspace(xmin - pad, xmax + pad, n_grid)
    gy = np.linspace(ymin - pad, ymax + pad, n_grid)
    GX, GY = np.meshgrid(gx, gy)

    # Interpolate velocity onto regular grid
    GU = griddata(embed, vel[:, 0], (GX, GY), method="linear")
    GV = griddata(embed, vel[:, 1], (GX, GY), method="linear")

    if stream:
        # Mask grid points with no data (outside convex hull)
        mask = np.isnan(GU) | np.isnan(GV)
        GU[mask] = 0.0; GV[mask] = 0.0
        speed = np.sqrt(GU**2 + GV**2)
        speed[mask] = 0.0

        ax.streamplot(gx, gy, GU, GV,
                      color=speed, cmap="Greys", linewidth=0.9,
                      density=stream_density, arrowsize=1.2,
                      norm=plt.Normalize(speed[~mask].min(), speed[~mask].max()))
    else:
        # Quiver fallback
        gx2, gy2, gu2, gv2 = _grid_quiver(embed, vel, n_bins=20, min_cells=min_cells)
        sc2 = np.percentile(np.sqrt(gu2**2 + gv2**2) + 1e-8, 70)
        ax.quiver(gx2, gy2, gu2 / sc2, gv2 / sc2,
                  alpha=0.85, color="black", scale=12, width=0.004, headwidth=4)

    xlabel = "PC 1" if basis == "X_pca" else "UMAP 1"
    ylabel = "PC 2" if basis == "X_pca" else "UMAP 2"
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    ax.set_title("Inferred drift field", fontsize=10)

    if standalone:
        plt.tight_layout()
        if save:
            plt.savefig(save, dpi=300, bbox_inches="tight")
        plt.show()
    return ax


def sensitivity(
    adata,
    key: str = "scjdo",
    ax: Optional[plt.Axes] = None,
    save: Optional[str] = None,
):
    """Plot max real eigenvalue (local sensitivity) across pseudotime."""
    t_np    = _get(adata, key, "t_centers")
    eig     = _get(adata, key, "max_real_eig")
    stable  = eig <= -0.05
    sens    = eig >= +0.05

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(6, 3))

    ax.axhline(0, color="gray", lw=0.8, ls="--")
    ax.axhline(+0.05, color="#E63946", lw=0.7, ls=":", alpha=0.8, label="Sensitive (≥+0.05)")
    ax.axhline(-0.05, color="#457B9D", lw=0.7, ls=":", alpha=0.8, label="Stable (≤-0.05)")
    ax.plot(t_np, eig, color="black", lw=1.8, label="Max Re(λ)")
    ax.fill_between(t_np, eig, 0, where=sens,   alpha=0.30, color="#E63946")
    ax.fill_between(t_np, eig, 0, where=stable,  alpha=0.20, color="#457B9D")
    ax.set_xlabel("Pseudotime"); ax.set_ylabel("Max Re(λ)")
    ax.set_title("Local sensitivity", fontsize=10)
    ax.legend(fontsize=8)

    if standalone:
        plt.tight_layout()
        if save:
            plt.savefig(save, dpi=300, bbox_inches="tight")
        plt.show()
    return ax


def archetypes(
    adata,
    key: str = "scjdo",
    ax: Optional[plt.Axes] = None,
    save: Optional[str] = None,
):
    """Plot archetype temporal activation profiles."""
    t_np      = _get(adata, key, "t_centers")
    act_norm  = _get(adata, key, "act_norm")
    K         = act_norm.shape[1]
    labels    = [f"A{k+1}" for k in range(K)]
    colors    = _ARCH_COLORS[:K]

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(6.5, 3.5))

    for k in range(K):
        ax.plot(t_np, act_norm[:, k], color=colors[k], lw=2, label=labels[k])
    ax.axhline(0, color="gray", lw=0.5, ls="--")
    ax.set_xlabel("Pseudotime"); ax.set_ylabel("Normalised activation")
    ax.set_title(f"Archetype activation profiles (K={K})", fontsize=10)
    ax.legend(ncol=min(K, 5), fontsize=9, loc="upper right")

    if standalone:
        plt.tight_layout()
        if save:
            plt.savefig(save, dpi=300, bbox_inches="tight")
        plt.show()
    return ax


def coordination(
    adata,
    key: str = "scjdo",
    axes: Optional[tuple] = None,
    save: Optional[str] = None,
):
    """Plot archetype coordination: correlation heatmap + activation timeline."""
    corr_mat = _get(adata, key, "corr_mat")
    t_np     = _get(adata, key, "t_centers")
    act_norm = _get(adata, key, "act_norm")
    K        = corr_mat.shape[0]
    labels   = [f"A{k+1}" for k in range(K)]
    colors   = _ARCH_COLORS[:K]

    standalone = axes is None
    if standalone:
        fig, (ax_h, ax_t) = plt.subplots(1, 2, figsize=(11, 4))
    else:
        ax_h, ax_t = axes

    # Heatmap
    im = ax_h.imshow(corr_mat, vmin=-1, vmax=1, cmap="RdBu_r", aspect="auto")
    ax_h.set_xticks(range(K)); ax_h.set_xticklabels(labels, fontsize=9)
    ax_h.set_yticks(range(K)); ax_h.set_yticklabels(labels, fontsize=9)
    for i in range(K):
        for j in range(K):
            c = corr_mat[i, j]
            ax_h.text(j, i, f"{c:+.2f}", ha="center", va="center", fontsize=8,
                      color="white" if abs(c) > 0.6 else "black")
    plt.colorbar(im, ax=ax_h, label="Temporal r", fraction=0.046, pad=0.04)
    ax_h.set_title("Archetype coordination", fontsize=10)

    # Timeline
    for k in range(K):
        ax_t.plot(t_np, act_norm[:, k], color=colors[k], lw=2, label=labels[k])
    ax_t.axhline(0, color="gray", lw=0.5, ls="--")
    # Mark sequential handoffs
    from itertools import combinations
    for i, j in combinations(range(K), 2):
        if corr_mat[i, j] < -0.5:
            diff  = act_norm[:, i] - act_norm[:, j]
            cross = np.where(np.diff(np.sign(diff)))[0]
            for c in cross:
                ax_t.axvline(t_np[c], color="#555", ls="--", lw=1, alpha=0.6)
    ax_t.set_xlabel("Pseudotime"); ax_t.set_ylabel("Activation")
    ax_t.set_title("Coordination timeline", fontsize=10)
    ax_t.legend(ncol=min(K, 5), fontsize=8)

    if standalone:
        plt.tight_layout()
        if save:
            plt.savefig(save, dpi=300, bbox_inches="tight")
        plt.show()
    return ax_h, ax_t


def instability_genes(
    adata,
    key: str = "scjdo",
    n_genes: int = 10,
    sensitivity_threshold: float = 0.05,
    per_archetype: bool = True,
    save: Optional[str] = None,
):
    """
    Plot and tabulate the top genes driving local instability across pseudotime.

    For each sensitive pseudotime window (Re(λ_max) > sensitivity_threshold),
    the eigenvector associated with the maximum real eigenvalue is projected to
    gene space via PCA loadings. Genes with the highest loadings on this
    "unstable direction" are those whose perturbation is most amplified by the
    drift field at that point in the trajectory.

    Parameters
    ----------
    n_genes              : Number of top genes to show in plots and table.
    sensitivity_threshold: Min Re(λ_max) to consider a window sensitive.
    per_archetype        : Also show top instability genes per archetype.

    Returns
    -------
    table : pandas.DataFrame
        Columns: gene, mean_instability_score, peak_pseudotime,
                 primary_archetype, [archetype activation at peak]
    """
    import pandas as pd

    res        = _get(adata, key, "instability_scores")   # (n_windows, n_genes)
    t_np       = _get(adata, key, "t_centers")
    max_eig    = _get(adata, key, "max_real_eig")
    act_norm   = _get(adata, key, "act_norm")
    gene_names = _get(adata, key, "gene_names")
    top_genes  = _get(adata, key, "top_instability_genes")
    top_scores = _get(adata, key, "top_instability_scores")
    arch_genes = _get(adata, key, "arch_instability_genes")

    if not gene_names:
        raise ValueError("No gene names found. Re-run sjd.tl.fit_drift on an adata with PCA loadings.")

    K          = act_norm.shape[1]
    arch_labels = [f"A{k+1}" for k in range(K)]
    sens_mask  = max_eig > sensitivity_threshold
    top_n      = top_genes[:n_genes]

    # ── Figure layout ──────────────────────────────────────────────────────
    n_rows = 3 if per_archetype else 2
    fig, axes = plt.subplots(n_rows, 1,
                             figsize=(10, 4 * n_rows),
                             gridspec_kw={"hspace": 0.45})

    # ── Panel 1: sensitivity curve + sensitive windows shaded ──────────────
    ax0 = axes[0]
    ax0.fill_between(t_np, 0, max_eig, where=sens_mask,
                     alpha=0.25, color="#E63946", label="Sensitive")
    ax0.fill_between(t_np, max_eig, 0, where=~sens_mask & (max_eig < 0),
                     alpha=0.15, color="#457B9D", label="Stable")
    ax0.plot(t_np, max_eig, color="black", lw=1.5)
    ax0.axhline(0, color="gray", lw=0.8, ls="--")
    ax0.axhline(sensitivity_threshold, color="#E63946", lw=0.8, ls=":", alpha=0.7)
    ax0.set_xlabel("Pseudotime"); ax0.set_ylabel("Max Re(λ)")
    ax0.set_title("Local sensitivity — sensitive windows (shaded red) drive instability", fontsize=10)
    ax0.legend(fontsize=8, loc="upper right")

    # ── Panel 2: top instability genes across pseudotime ───────────────────
    ax1    = axes[1]
    cmap2  = plt.get_cmap("tab10", n_genes)
    g_idx  = [gene_names.index(g) for g in top_n if g in gene_names]

    for rank, gi in enumerate(g_idx):
        scores = res[:, gi].copy()
        # Zero out non-sensitive windows
        scores[~sens_mask] = np.nan
        ax1.plot(t_np, scores, color=cmap2(rank), lw=1.8,
                 label=gene_names[gi], alpha=0.85)

    # Light sensitivity background
    ax1.fill_between(t_np, ax1.get_ylim()[0] if ax1.get_ylim()[0] < 0 else -0.05,
                     0.05, where=sens_mask, alpha=0.06, color="#E63946")
    ax1.axhline(0, color="gray", lw=0.5, ls="--")
    ax1.set_xlabel("Pseudotime")
    ax1.set_ylabel("Instability gene score")
    ax1.set_title(f"Top {n_genes} instability-driving genes across pseudotime\n"
                  "(score = projection onto max-eigenvalue eigenvector; shown only in sensitive windows)",
                  fontsize=10)
    ax1.legend(ncol=min(n_genes, 5), fontsize=7, loc="upper right")

    # ── Panel 3: heatmap of top genes per archetype ────────────────────────
    if per_archetype:
        ax2 = axes[2]
        # Build matrix: top_genes × archetypes, value = mean instability score
        heatmap_data = np.zeros((n_genes, K), dtype=np.float32)
        row_labels   = []
        for row, gene in enumerate(top_n):
            if gene not in gene_names:
                continue
            row_labels.append(gene)
            gi = gene_names.index(gene)
            for k in range(K):
                thresh    = np.quantile(act_norm[:, k], 0.75)
                arch_mask = (act_norm[:, k] >= thresh) & sens_mask
                if arch_mask.sum() > 0:
                    heatmap_data[row, k] = float(np.abs(res[arch_mask, gi]).mean())

        im2 = ax2.imshow(heatmap_data[:len(row_labels)], aspect="auto",
                         cmap="YlOrRd", interpolation="nearest")
        ax2.set_xticks(range(K)); ax2.set_xticklabels(arch_labels, fontsize=9)
        ax2.set_yticks(range(len(row_labels))); ax2.set_yticklabels(row_labels, fontsize=8)
        plt.colorbar(im2, ax=ax2, label="Mean |instability score|", fraction=0.03, pad=0.02)
        ax2.set_xlabel("Archetype")
        ax2.set_title("Mean instability score per gene × archetype\n"
                      "(high = gene drives instability when this archetype is active)", fontsize=10)

    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=300, bbox_inches="tight")
    plt.show()

    # ── Build summary table ────────────────────────────────────────────────
    rows = []
    for rank, (gene, score) in enumerate(zip(top_genes[:n_genes],
                                              top_scores[:n_genes])):
        if gene not in gene_names:
            continue
        gi = gene_names.index(gene)

        # Pseudotime of peak instability score
        g_scores_sens       = res[:, gi].copy()
        g_scores_sens[~sens_mask] = 0
        peak_pt             = float(t_np[np.argmax(np.abs(g_scores_sens))])

        # Primary archetype: which archetype is most active at peak window?
        peak_win            = int(np.argmax(np.abs(g_scores_sens)))
        primary_arch        = int(np.argmax(act_norm[peak_win])) + 1

        rows.append({
            "rank":                   rank + 1,
            "gene":                   gene,
            "mean_instability_score": round(float(score), 4),
            "peak_pseudotime":        round(peak_pt, 3),
            "primary_archetype":      f"A{primary_arch}",
        })

    table = pd.DataFrame(rows)
    print(table.to_string(index=False))
    return table


def summary_figure(
    adata,
    key: str = "scjdo",
    basis: str = "X_pca",
    save: Optional[str] = None,
):
    """
    Assemble the four-panel Figure 3 layout:
    (a) drift field, (b) local sensitivity, (c) archetypes, (d) coordination.
    """
    fig = plt.figure(figsize=(13, 10))
    gs  = gridspec.GridSpec(2, 2, hspace=0.40, wspace=0.35,
                            left=0.07, right=0.97, top=0.92, bottom=0.08)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    drift_field(adata, key=key, basis=basis, ax=ax_a)
    sensitivity(adata, key=key, ax=ax_b)
    archetypes(adata, key=key, ax=ax_c)
    # Coordination heatmap only in panel d
    corr_mat = _get(adata, key, "corr_mat")
    K        = corr_mat.shape[0]
    labels   = [f"A{k+1}" for k in range(K)]
    im = ax_d.imshow(corr_mat, vmin=-1, vmax=1, cmap="RdBu_r", aspect="auto")
    ax_d.set_xticks(range(K)); ax_d.set_xticklabels(labels, fontsize=9)
    ax_d.set_yticks(range(K)); ax_d.set_yticklabels(labels, fontsize=9)
    for i in range(K):
        for j in range(K):
            c = corr_mat[i, j]
            ax_d.text(j, i, f"{c:+.2f}", ha="center", va="center", fontsize=8,
                      color="white" if abs(c) > 0.6 else "black")
    plt.colorbar(im, ax=ax_d, label="Temporal r", fraction=0.046, pad=0.04)
    ax_d.set_title("d  |  Archetype coordination", fontweight="bold", loc="left", fontsize=10)

    for ax, label in zip([ax_a, ax_b, ax_c], ["a", "b", "c"]):
        ax.set_title(f"{label}  |  " + ax.get_title(), fontweight="bold", loc="left", fontsize=10)

    fig.suptitle("scJDO — operator-level view of single-cell dynamics",
                 fontsize=11, fontweight="bold")

    if save:
        fig.savefig(save, dpi=300, bbox_inches="tight")
        print(f"Saved: {save}")
    plt.show()
    return fig
