# Functional Pathway Analysis of scJDO / scOpAtlas Instability Drivers

This guide shows how to take the **instability-driver genes** extracted from a fitted scJDO model and analyze them with either over-representation analysis (Enrichr) or ranked-list enrichment analysis (GSEA). It is designed for the unified scJDO + scOpAtlas workflow.

> **Interpretation first.** A scJDO instability-driver score is not automatically a differential-expression effect size. It measures association with an unstable direction in the learned drift field. Enrichment analysis therefore asks which annotated functions are concentrated among genes associated with an inferred dynamical transition; it does not by itself prove that a pathway caused that transition.

## 1. Choose the correct enrichment question

| Question | Recommended analysis | Required input | Interpretation |
|---|---|---|---|
| Which functions are over-represented among a short, high-confidence driver list? | Enrichr / over-representation analysis | One gene list plus a measured-gene background | Functions concentrated among selected instability drivers |
| Are pathway members preferentially concentrated across the complete instability ranking? | Pre-ranked GSEA | A complete, ordered gene-score vector | Functions distributed toward the high-driver end of the score distribution |
| Do different dynamic archetypes implicate distinct programs? | Run either approach **separately for each archetype** | One list/ranking per archetype | Archetype-specific dynamical programs |
| Does a pathway distinguish a condition or perturbation? | Differential-expression or differential-abundance analysis, optionally alongside scJDO | Condition-aware design | A condition effect; do not infer it from driver ranking alone |

Enrichr accepts an input list of gene symbols and tests it against curated libraries; it also supports an explicit background through its interface/API. [1] GSEA PreRank evaluates whether predefined gene sets are concentrated at the top or bottom of a user-supplied ranked list. [2] GSEApy provides Python interfaces for Enrichr, over-representation analysis, and pre-ranked GSEA. [3]

## 2. Export scJDO instability drivers

After fitting the drift model, extract archetype-level drivers. This returns a table with a gene, rank, instability score, archetype, and sensitivity-related summary fields.

```python
import scjdo as sjd

# Requires adata.uns['scjdo'] created by sjd.tl.fit_drift(...).
driver_table = sjd.tl.get_instability_genes(
    adata,
    key="scjdo",
    n_genes=100,
    activation_threshold=0.5,
    sensitivity_threshold=0.05,
    min_sensitive_fraction=0.10,
)

driver_table.to_csv("results/instability_driver_genes_by_archetype.csv", index=False)
print(driver_table.head())
```

A missing or empty table is informative: it means no archetype met the selected instability/sensitivity criteria. Do **not** lower thresholds merely to force a gene list. Instead, inspect the drift-quality report, training convergence, and whether the dataset contains a transition supported by the data.

## 3. Quality control before enrichment

Perform the following checks before sending genes to any pathway tool.

| Check | Why it matters | Practical action |
|---|---|---|
| Drift quality gate | Poorly supported dynamics can yield spurious driver rankings. | Inspect `adata.uns['drift_quality']`; do not use forced atlas results as confirmatory evidence. |
| Gene identifiers | Enrichment libraries require matching identifiers. | Use official gene symbols consistently; remove null, duplicate, and ambiguous aliases. |
| Organism | Human and mouse collections differ. | Select organism-matched Enrichr/MSigDB collections or map orthologs explicitly. |
| Background | ORA significance is biased if unmeasured genes are treated as testable. | Use genes measured and eligible in the assay, not the whole genome by default. |
| Minimum list size | Very short lists are unstable; very long lists become non-specific. | Report list size and sensitivity threshold with every result. |
| Multiple archetypes | A pooled list can hide lineage-specific drivers. | Analyze each archetype separately, then compare pathways. |

Construct an assay-aware background from the genes that were tested by the scJDO analysis. For a PCA-based workflow, this is normally the post-HVG feature set stored in `adata.var_names`.

```python
import pandas as pd

background = (
    pd.Index(adata.var_names)
      .astype(str)
      .str.strip()
      .drop_duplicates()
      .tolist()
)
```

## 4. Enrichr over-representation analysis

Use Enrichr when the scientific question is about the top selected genes. GSEApy can submit a list and preserve the results in a local table. Its documentation shows that `gp.enrichr()` accepts gene lists, library names, organism choice, and a custom background. [3]

### 4.1 Run an archetype-specific Enrichr analysis

```python
import pandas as pd
import gseapy as gp

# Confirm currently available libraries rather than assuming a fixed release name.
mouse_libraries = gp.get_library_name(organism="Mouse")
print(mouse_libraries[:20])

# Pick library names present in the current service response.
selected_libraries = [
    name for name in [
        "GO_Biological_Process_2023",
        "KEGG_2021_Mouse",
        "WikiPathway_2024_Mouse",
    ]
    if name in mouse_libraries
]
if not selected_libraries:
    raise RuntimeError("No requested mouse Enrichr libraries found; inspect mouse_libraries.")

for archetype, subset in driver_table.groupby("archetype", observed=True):
    genes = (
        subset.sort_values("instability_score", ascending=False)["gene"]
              .astype(str).str.strip()
              .drop_duplicates()
              .tolist()
    )

    if len(genes) < 5:
        print(f"Skipping {archetype}: only {len(genes)} unique driver genes.")
        continue

    result = gp.enrichr(
        gene_list=genes,
        gene_sets=selected_libraries,
        organism="Mouse",
        background=background,
        outdir=None,
    )

    table = result.results.sort_values("Adjusted P-value")
    table.to_csv(f"results/enrichr_{archetype}.csv", index=False)
    print(f"\n{archetype}: {len(genes)} genes")
    print(table[["Gene_set", "Term", "Adjusted P-value", "Odds Ratio", "Genes"]].head(10))
```

### 4.2 Report Enrichr results responsibly

Report the gene-list construction rule, organism, libraries, explicit background, number of input genes, library version/date, multiple-testing statistic, and genes overlapping each term. Use adjusted p-values as the primary significance measure and treat overlapping pathway labels as correlated annotations rather than independent discoveries.

The Enrichr interface explicitly supports input gene lists and an optional background; supply the background whenever the assay measured only a subset of genes. [1]

## 5. Pre-ranked GSEA on the complete instability score

GSEA is preferable when a complete ranking exists, because it avoids an arbitrary top-*N* cutoff. GSEA PreRank expects a ranked input containing unique identifiers and unique ranking values; duplicated scores can produce arbitrary ordering. [2]

### 5.1 Build a complete instability ranking

The code below derives a global ranking from the mean driver score across sensitive windows. It is an example; for an archetype-specific GSEA, replace the aggregation with the matching archetype score vector if that full vector was retained.

```python
import numpy as np
import pandas as pd

result = adata.uns["scjdo"]
window_scores = np.asarray(result["instability_scores"], dtype=float)
gene_names = pd.Index(result["gene_names"]).astype(str)
max_real_eig = np.asarray(result["max_real_eig"], dtype=float)

sensitivity_threshold = 0.05
sensitive = max_real_eig > sensitivity_threshold
if not sensitive.any():
    raise RuntimeError("No sensitive windows; do not run driver-score GSEA.")

# Higher score means stronger association with the selected unstable windows.
ranked = pd.Series(window_scores[sensitive].mean(axis=0), index=gene_names)
ranked = ranked.groupby(level=0).max().dropna().sort_values(ascending=False)

# GSEA PreRank requires unique ranking values. A deterministic tiny tie-breaker
# preserves rank ordering while avoiding arbitrary duplicate-score ordering.
ranked = ranked + np.linspace(0.0, 1e-12, len(ranked), endpoint=False)
ranked.index.name = "gene"
ranked.name = "instability_score"
ranked.to_csv("results/global_instability_driver.rnk", sep="\t", header=False)
```

### 5.2 Run pre-ranked GSEA with GSEApy

```python
import gseapy as gp

# For mouse data, use a mouse MSigDB collection or a custom mouse GMT.
# The GSEApy MSigDB API can fetch an organism-specific GMT collection.
msig = gp.Msigdb()
mouse_hallmark = msig.get_gmt(category="mh.all", dbver="2023.1.Mm")

pre_res = gp.prerank(
    rnk=ranked,
    gene_sets=mouse_hallmark,
    min_size=15,
    max_size=500,
    permutation_num=1_000,
    seed=17,
    outdir="results/gsea_global_instability",
    verbose=True,
)

columns = ["Term", "ES", "NES", "NOM p-val", "FDR q-val", "Lead_genes"]
summary = pre_res.res2d.sort_values("FDR q-val")
summary.to_csv("results/gsea_global_instability_summary.csv", index=False)
print(summary[columns].head(15))
```

GSEApy documents `prerank` as a single-ranked-list analysis and supports gene-set input as a collection name, local GMT, or in-memory dictionary. [3] The GSEA PreRank specification recommends unique gene identifiers and unique ranking values; it uses gene-set permutations and commonly uses 1,000 permutations for a completed analysis. [2]

### 5.3 Interpret the direction correctly

For a conventional signed differential-expression statistic, positive and negative normalized enrichment scores have a direct up/down interpretation. scJDO instability scores are often non-negative magnitudes. In that case, a high positive enrichment score means that a pathway’s genes are concentrated among stronger **instability-associated** genes; a low-end enrichment signal means those genes are relatively depleted from that ranking. It does **not** mean the pathway is transcriptionally down-regulated.

To obtain a signed biological interpretation, construct a score that retains a justified direction—for example, a signed projection of genes onto a specified unstable mode—and validate that direction independently with expression, perturbation, or temporal evidence.

## 6. Compare enrichment across operator regimes or archetypes

A practical analysis matrix is:

| Analysis unit | Driver input | Recommended test | Primary comparison |
|---|---|---|---|
| One unstable archetype | Top drivers + full score ranking | Enrichr + pre-ranked GSEA | Which functions define this instability mode? |
| Multiple archetypes | Separate lists/rankings | Per-archetype enrichment | Which pathways are specific versus shared? |
| Operator regime | Differential expression or differential score conditioned on covariates | DE plus GSEA | Which pathways differ between stable/plastic/unstable states? |
| Experimental condition | Condition-aware model within an operator regime | DE plus GSEA | Does perturbation alter a regime’s program? |

Use the same background and gene-set collection for a comparison across archetypes. Summarize shared pathways with leading-edge gene overlap, not only term-name overlap. A robust finding should recur across model seeds, relevant samples, and reasonable threshold choices.

## 7. Recommended validation ladder

1. **Computational robustness:** re-fit scJDO with different seeds/subsamples and test whether driver rankings and enriched pathways recur.
2. **Representation robustness:** compare PCA versus a batch-aware latent representation when batch effects are plausible.
3. **Biological coherence:** compare against known lineage, metabolic, or perturbation programs without using those programs to define the driver score.
4. **Independent evidence:** test key leading-edge genes with differential expression, protein, spatial, perturbation, or external dataset evidence.
5. **Causal follow-up:** perturb candidate drivers only after the preceding evidence supports a stable, context-specific hypothesis.

## 8. Common failure modes

| Failure mode | Why it misleads | Correct response |
|---|---|---|
| Enriching a forced/weak driver list | The ranking may reflect model noise, not a real transition. | Stop and improve/validate the drift model. |
| Using a whole-genome background for an HVG-filtered assay | Inflates apparent over-representation. | Use measured/tested genes as background. |
| Mixing human libraries with mouse symbols | Produces poor overlap or invalid interpretation. | Use organism-matched libraries or explicit ortholog mapping. |
| Treating pathway labels as independent discoveries | Libraries and terms overlap extensively. | Collapse redundant terms and inspect leading-edge genes. |
| Calling a driver score “up-regulation” | Instability is a dynamical association, not necessarily expression change. | Use neutral language unless signed expression evidence is available. |
| Pooling all archetypes before analysis | Obscures mode-specific biology. | Analyze archetypes separately, then compare. |

## References

[1] [Enrichr: Ma’ayan Lab](https://maayanlab.cloud/enrichr/)

[2] [GSEAPreranked: GenePattern documentation](https://www.genepattern.org/modules/docs/GSEAPreranked/1/)

[3] [GSEApy tutorial and API examples](https://gseapy.readthedocs.io/en/latest/gseapy_example.html)

[4] [Subramanian et al., 2005: Gene set enrichment analysis](https://doi.org/10.1073/pnas.0506580102)
