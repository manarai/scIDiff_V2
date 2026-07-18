"""Rerun the multiome pipeline WITHOUT per-modality CC regression.

Purpose. The manuscript's Supplementary Note claims that not regressing
S/G2M scores from RNA + ATAC before integration leaves the ExcNeuron top-9
dominated by Cdc20 / Ccnb1 / Cenpf. The shipped `adata_multiome_fa.h5ad`
was built with `DO_CC_REGRESSION = True` and cannot be used to test that
claim (all three FA latents in the input are already orthogonal to S/G2M).

This script reruns the essentials of `Figure6_02_multiome_FA_integration.ipynb`
with `DO_CC_REGRESSION = False`, then borrows the ExcNeuron branch
definition and Palantir pseudotime from the corrected pipeline (for a fair
apples-to-apples comparison — same cells, same drift-training config,
different latent), fits scJDO drift, and reports the ExcNeuron top-9.

Outputs
  analysis/results/supp_multiome/uncorrected_pipeline/
    - adata_multiome_fa_uncorrected.h5ad
    - instability_genes_ExcNeuron_uncorrected.csv
    - top9_comparison.json  (uncorrected top-9 vs currently shipped
                             pipelines' top-9 + mitotic-gene overlap flags)
"""
from __future__ import annotations

import json
import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
from sklearn.decomposition import FactorAnalysis
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, "/Users/terooatt/Downloads/scJDO")
import scjdo as sjd

warnings.filterwarnings("ignore")

# ── Paths ────────────────────────────────────────────────────────────────
DATA_DIR = "/Users/terooatt/Documents/Project_scQDiff/02_scQDiff/scIDIFF_anndata/data/multiomic/filtered_feature_bc_matrix"
FIG6_DIR = "/Users/terooatt/Downloads/scJDO/Figures_notebook/results/figure6_multiome_fa"
DRIFT_CC_H5AD = "/Users/terooatt/Downloads/scJDO/Figures_notebook/results/figure6_multiome_fa_drift_cc/adata_multiome_drift_cc.h5ad"
OUTDIR = Path("/Users/terooatt/Downloads/scJDO/analysis/results/supp_multiome/uncorrected_pipeline")
OUTDIR.mkdir(parents=True, exist_ok=True)

# ── Config (mirrors Figure6_02) ──────────────────────────────────────────
N_HVG = 2000
N_VAR_PEAKS = 10000
N_LATENT = 30
FA_MAX_ITER = 1000
TFIDF_SCALE = 1e4
MIN_CELLS_PEAK = 25
MIN_CELLS_GENE = 10
SEED = 42

N_ARCHETYPES = 5
N_EPOCHS = 5000
BIAS_STRENGTH = 1.5


def tfidf_transform(X, scale=TFIDF_SCALE):
    X = X.tocsr().astype(np.float32)
    n_cells = X.shape[0]
    cell_totals = np.asarray(X.sum(axis=1)).ravel()
    cell_totals[cell_totals == 0] = 1.0
    peak_df = np.asarray((X > 0).sum(axis=0)).ravel()
    idf = np.log(1.0 + n_cells / (1.0 + peak_df)).astype(np.float32)
    inv_totals = sp.diags(1.0 / cell_totals)
    tf = inv_totals @ X
    tfidf = tf.multiply(idf[None, :])
    tfidf.data = np.log1p(tfidf.data * scale)
    return tfidf.tocsr(), idf


print("=" * 72)
print("UNCORRECTED multiome pipeline — no S/G2M regression at any stage")
print("=" * 72)

# ── Load joint 10x matrix ────────────────────────────────────────────────
t0 = time.time()
ad_full = sc.read_10x_mtx(DATA_DIR, gex_only=False, cache=True)
ad_full.var_names_make_unique()
print(f"Joint 10x: {ad_full.shape}  ({time.time()-t0:.1f}s)")
print("Feature types:")
print(ad_full.var["feature_types"].value_counts().to_string())

# Split into RNA + ATAC
rna_mask = (ad_full.var["feature_types"] == "Gene Expression").values
atac_mask = (ad_full.var["feature_types"] == "Peaks").values
ad_rna = ad_full[:, rna_mask].copy()
ad_atac = ad_full[:, atac_mask].copy()
print(f"RNA:  {ad_rna.shape}")
print(f"ATAC: {ad_atac.shape}")

# ── Step 2: RNA preprocessing ────────────────────────────────────────────
sc.pp.filter_genes(ad_rna, min_cells=MIN_CELLS_GENE)
ad_rna.layers["raw_counts"] = ad_rna.X.copy()
sc.pp.normalize_total(ad_rna, target_sum=1e4)
sc.pp.log1p(ad_rna)
sc.pp.highly_variable_genes(ad_rna, n_top_genes=N_HVG, flavor="seurat")

rna_hvg_mask = ad_rna.var["highly_variable"].values
X_rna = ad_rna[:, rna_hvg_mask].X
X_rna = (X_rna.toarray() if sp.issparse(X_rna) else np.asarray(X_rna)).astype(np.float32)
print(f"RNA HVG: {X_rna.shape}")

# ── Step 3: ATAC TF-IDF + variable peaks ─────────────────────────────────
sc.pp.filter_genes(ad_atac, min_cells=MIN_CELLS_PEAK)
ad_atac.layers["raw_counts"] = ad_atac.X.copy()
X_atac_full, _idf = tfidf_transform(ad_atac.X)
mean_sq = np.asarray(X_atac_full.multiply(X_atac_full).mean(axis=0)).ravel()
mean = np.asarray(X_atac_full.mean(axis=0)).ravel()
peak_var = mean_sq - mean**2
top_idx = np.argsort(peak_var)[::-1][:N_VAR_PEAKS]
atac_var_mask = np.zeros(ad_atac.n_vars, dtype=bool)
atac_var_mask[top_idx] = True
X_atac = X_atac_full[:, atac_var_mask].toarray().astype(np.float32)
print(f"ATAC variable peaks: {X_atac.shape}")

# ── Step 3b: SKIP CC regression ──────────────────────────────────────────
print("\n[SKIP] Step 3b — DO_CC_REGRESSION = False (this is the whole point)")

# ── Step 4: balance + concat ─────────────────────────────────────────────
scaler_rna = StandardScaler(with_mean=True, with_std=True).fit(X_rna)
scaler_atac = StandardScaler(with_mean=True, with_std=True).fit(X_atac)
X_rna_s = scaler_rna.transform(X_rna).astype(np.float32)
X_atac_s = scaler_atac.transform(X_atac).astype(np.float32)
w_rna = 1.0 / np.sqrt(X_rna_s.shape[1])
w_atac = 1.0 / np.sqrt(X_atac_s.shape[1])
X_joint = np.concatenate([X_rna_s * w_rna, X_atac_s * w_atac], axis=1)
n_rna_feats = X_rna_s.shape[1]
print(f"Joint matrix: {X_joint.shape}")

# ── Step 5: joint FA (uncorrected) ───────────────────────────────────────
t0 = time.time()
fa_joint = FactorAnalysis(n_components=N_LATENT, random_state=SEED,
                          max_iter=FA_MAX_ITER)
Z_joint = fa_joint.fit_transform(X_joint).astype(np.float32)
print(f"Joint FA (uncorrected) fit in {time.time()-t0:.1f}s")

# Sanity check: does the UNCORRECTED latent carry CC signal?
sc.tl.score_genes_cell_cycle(
    ad_rna,
    s_genes=[g for g in [
        "Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6",
        "Cdca7","Dtl","Prim1","Uhrf1","Mlf1ip","Hells","Rfc2","Rpa2","Nasp",
        "Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2",
        "Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2",
        "Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8",
    ] if g in ad_rna.var_names],
    g2m_genes=[g for g in [
        "Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2",
        "Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Fam64a","Smc4","Ccnb2",
        "Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1",
        "Kif20b","Hjurp","Cdca3","Hn1","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1",
        "Ncapd2","Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1",
        "Anln","Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa",
    ] if g in ad_rna.var_names],
)
from scipy.stats import pearsonr
print("\nUncorrected latent — |Pearson| with S/G2M (first 5 factors):")
for k in range(5):
    r_s = pearsonr(Z_joint[:, k], ad_rna.obs["S_score"].values)[0]
    r_g = pearsonr(Z_joint[:, k], ad_rna.obs["G2M_score"].values)[0]
    print(f"  factor {k}: |r_S|={abs(r_s):.3f}, |r_G2M|={abs(r_g):.3f}")

# ── Borrow branch definitions + pseudotime from the shipped CC pipeline ─
print("\nLoading shipped CC-corrected h5ad for branch labels + pseudotime …")
ad_cc = sc.read_h5ad(DRIFT_CC_H5AD)
common = ad_rna.obs_names.intersection(ad_cc.obs_names)
print(f"  cells in common: {len(common)}")

ad_joint = ad_rna.copy()
ad_joint.obsm["X_fa_joint"] = Z_joint

# Zero-padded RNA loadings for scJDO
rna_loadings = fa_joint.components_[:, :n_rna_feats].T.astype(np.float32)
full_rna_loadings = np.zeros((ad_rna.n_vars, N_LATENT), dtype=np.float32)
full_rna_loadings[rna_hvg_mask] = rna_loadings
ad_joint.varm["PCs"] = full_rna_loadings

# Transfer branch masks + pseudotime + fate labels from CC h5ad
for col in ["branch_masks_cc"]:
    if col in ad_cc.obsm:
        m_cc = ad_cc.obsm[col].loc[common]
        m_new = pd.DataFrame(False, index=ad_joint.obs_names,
                             columns=m_cc.columns)
        m_new.loc[common] = m_cc.values
        ad_joint.obsm[col] = m_new
for col in ["palantir_pseudotime_cc", "cell_fate_cc"]:
    if col in ad_cc.obs:
        vals = pd.Series(np.nan, index=ad_joint.obs_names, dtype=object)
        vals.loc[common] = ad_cc.obs[col].loc[common].values
        if col == "palantir_pseudotime_cc":
            ad_joint.obs[col] = vals.astype(float)
        else:
            ad_joint.obs[col] = pd.Categorical(vals)

# Drop cells not in the CC pipeline (they don't have labels/pseudotime)
mask = ad_joint.obs_names.isin(common)
ad_joint = ad_joint[mask].copy()
print(f"  kept {ad_joint.n_obs} cells (matched to CC h5ad)")

# Reduce branch mask to same cells
for col in ["branch_masks_cc"]:
    if col in ad_joint.obsm:
        ad_joint.obsm[col] = ad_joint.obsm[col].loc[ad_joint.obs_names]

# ── Fit scJDO drift on the UNCORRECTED latent (X_fa_joint) ───────────────
BRANCH_NAMES = list(ad_joint.obsm["branch_masks_cc"].columns)
print(f"\nBranches (from CC pipeline, applied to UNCORRECTED latent): {BRANCH_NAMES}")

t0 = time.time()
branch_models = sjd.tl.fit_drift_branches(
    ad_joint,
    rep="X_fa_joint",
    branch_key="branch_masks_cc",
    branch_names=BRANCH_NAMES,
    time_key="palantir_pseudotime_cc",
    groupby="cell_fate_cc",
    progenitor_cluster="Progenitor",
    terminal_clusters={b: b for b in BRANCH_NAMES},
    bias_strength=BIAS_STRENGTH,
    n_archetypes=N_ARCHETYPES,
    n_epochs=N_EPOCHS,
    vel_scale=2.0,
    n_eff_min=50.0,
    seed=SEED,
    key_prefix="scjdo_uncorr",
)
print(f"Drift training (uncorrected): {time.time()-t0:.1f}s")

# ── Instability genes for each branch ────────────────────────────────────
print("\nComputing instability genes on the UNCORRECTED latent …")
instab_out = {}
for name in BRANCH_NAMES:
    key = f"scjdo_uncorr_{name}"
    if key not in ad_joint.uns:
        print(f"  {name}: {key} not in ad.uns — skipping")
        continue
    cell_idx = np.array(ad_joint.uns[key]["branch_cells"])
    ad_b = ad_joint[cell_idx].copy()
    ad_b.uns[key] = ad_joint.uns[key]
    try:
        table = sjd.pl.instability_genes(
            ad_b, key=key, n_genes=20,
            save=str(OUTDIR / f"instab_{name}_uncorrected.pdf"),
        )
    except Exception as e:
        print(f"  {name}: instability_genes failed ({type(e).__name__}: {e})")
        continue
    outpath = OUTDIR / f"instability_genes_{name}_uncorrected.csv"
    table.to_csv(outpath, index=False)
    gene_col = "gene" if "gene" in table.columns else table.columns[0]
    instab_out[name] = table[gene_col].tolist()
    print(f"  {name}: saved {outpath.name}")
    print(f"    top-9: {instab_out[name][:9]}")

# ── Comparison summary ───────────────────────────────────────────────────
mitotic = {"Cdc20", "Ccnb1", "Cenpf", "Mki67", "Top2a", "Aurkb",
           "Ccna2", "Ccnb2", "Kif20a", "Kif11", "Plk1", "Bub1", "Aurka"}

def load_top9(path):
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    return df.head(9)["gene"].tolist() if "gene" in df.columns else None

shipped_uncorrected_top9 = load_top9("/Users/terooatt/Downloads/scJDO/Figures_notebook/results/figure6_multiome_fa_drift/instability_genes_ExcNeuron.csv")
shipped_cc_top9 = load_top9("/Users/terooatt/Downloads/scJDO/Figures_notebook/results/figure6_multiome_fa_drift_cc/instability_genes_cc_ExcNeuron.csv")

def mitotic_hits(genes):
    return [g for g in (genes or []) if g in mitotic]

payload = {
    "purpose": "Test manuscript v38 Supplementary-Note claim that a truly uncorrected multiome pipeline (no per-modality CC regression) leaves ExcNeuron top-9 dominated by Cdc20/Ccnb1/Cenpf.",
    "run_config": {
        "DO_CC_REGRESSION": False,
        "n_cells_processed": int(ad_joint.n_obs),
        "n_HVG_rna": int(rna_hvg_mask.sum()),
        "n_var_peaks_atac": int(atac_var_mask.sum()),
        "n_latent": N_LATENT,
        "drift_n_epochs": N_EPOCHS,
    },
    "top9_uncorrected_this_run": {b: g[:9] for b, g in instab_out.items()},
    "top9_shipped_uncorrected_pipeline_ExcNeuron": shipped_uncorrected_top9,
    "top9_shipped_cc_pipeline_ExcNeuron": shipped_cc_top9,
    "mitotic_reference_set": sorted(mitotic),
    "mitotic_hits_ExcNeuron_uncorrected_this_run": mitotic_hits(instab_out.get("ExcNeuron", [])[:9]),
    "mitotic_hits_shipped_uncorrected_ExcNeuron": mitotic_hits(shipped_uncorrected_top9),
    "mitotic_hits_shipped_cc_ExcNeuron": mitotic_hits(shipped_cc_top9),
}

# Save adata for reference
try:
    ad_joint.write_h5ad(OUTDIR / "adata_multiome_fa_uncorrected.h5ad",
                        compression="gzip")
    print(f"\nSaved uncorrected h5ad: {OUTDIR / 'adata_multiome_fa_uncorrected.h5ad'}")
except Exception as e:
    print(f"h5ad save failed: {type(e).__name__}: {e}")

with open(OUTDIR / "top9_comparison.json", "w") as f:
    json.dump(payload, f, indent=2)
print(f"Saved comparison: {OUTDIR / 'top9_comparison.json'}")

print()
print("=" * 72)
print("VERDICT")
print("=" * 72)
if payload["mitotic_hits_ExcNeuron_uncorrected_this_run"]:
    print(f"  ExcNeuron top-9 (uncorrected this run) mitotic hits: "
          f"{payload['mitotic_hits_ExcNeuron_uncorrected_this_run']}")
    print("  → Manuscript claim supported: the uncorrected pipeline does show "
          "mitotic dominance in the ExcNeuron top-9.")
else:
    print("  ExcNeuron top-9 (uncorrected this run): no Cdc20 / Ccnb1 / Cenpf "
          "/ Mki67 / Top2a etc. hits.")
    print("  → Manuscript claim NOT reproducible even with the truly uncorrected "
          "pipeline (this dataset does not exhibit the described artifact).")
