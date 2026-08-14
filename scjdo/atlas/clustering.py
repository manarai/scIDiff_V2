"""
Operator-Based Clustering Utilities

Clustering methods that leverage operator space for improved cell state definition.

This module demonstrates that operator-based clustering provides information
beyond expression-based clustering, which is a key claim of the SCOPAtlas paper.
"""

import numpy as np
import scanpy as sc
import anndata as ad
from sklearn.metrics import adjusted_rand_score, silhouette_score, normalized_mutual_info_score
from sklearn.preprocessing import StandardScaler
from typing import Dict, List, Optional, Tuple
import warnings


class OperatorClustering:
    """
    Clustering utilities for operator space.
    
    This class provides methods for:
    1. Clustering cells using only operator features
    2. Joint clustering (expression + operator)
    3. Weighted combination of expression and operator spaces
    4. Clustering quality comparison
    
    Args:
        adata: AnnData object with operator metrics stored in .obs
        
    Example:
        >>> clusterer = OperatorClustering(adata)
        >>> clusterer.prepare_operator_features()
        >>> clusterer.cluster_operator_space()
        >>> clusterer.cluster_joint_space(alpha=0.5)
        >>> results = clusterer.compare_clustering_quality(
        ...     methods={'expression': 'leiden', 'operator': 'operator_clusters'},
        ...     celltype_key='cell_type'
        ... )
    """
    
    def __init__(self, adata: ad.AnnData):
        self.adata = adata
    
    def prepare_operator_features(
        self,
        metrics_keys: Optional[List[str]] = None,
        standardize: bool = True,
        obsm_key: str = 'X_operator',
    ) -> np.ndarray:
        """
        Prepare operator features matrix from metrics.

        Default features (Task 2.2): ['lambda_max_plus', 'stability_depth_q',
        'plasticity']. Two axes are dropped from the old four-tuple because
        they were redundant with the rest and z-scored Leiden treated them as
        independent, double-weighting the redundant axis:

          - `lambda_min_minus` (min-order statistic) is anti-correlated with
            `lambda_max_plus` by construction; use `stability_depth_q`, which
            is the dimension-stable quantile stability metric.
          - `stable_dim` (# of Re(λ) < 0) overlaps with `plasticity` for cells
            near neutrality (e.g. all-eigenvalues = -0.05 with ε = 0.1 scores
            both maximally plastic AND maximally buffered).

        Pass metrics_keys= to use a different feature set. To use spectrum-PC
        features instead, build them via OperatorEmbedding(method='spectrum')
        directly (requires store_eigenvalues=True on the atlas build).

        Args:
            metrics_keys: List of metric keys in adata.obs. If None, uses the
                three-feature default above and silently drops any keys not
                present (backwards-compat with older adata that lack
                stability_depth_q).
            standardize: Whether to standardize features
            obsm_key: Key to store operator features in adata.obsm

        Returns:
            operator_features: (n_cells, n_features) operator feature matrix

        Raises:
            ValueError: If any explicitly requested metric key is not found.
        """
        if metrics_keys is None:
            preferred = ['lambda_max_plus', 'stability_depth_q', 'plasticity']
            metrics_keys = [k for k in preferred if k in self.adata.obs]
            if not metrics_keys:
                raise ValueError(
                    "No default operator features found in adata.obs. "
                    "Run atlas.build() first."
                )
        features = []
        for key in metrics_keys:
            if key not in self.adata.obs:
                raise ValueError(f"Metric '{key}' not found in adata.obs. "
                               f"Run atlas.build() first to compute operator metrics.")
            features.append(self.adata.obs[key].values)
        
        operator_features = np.column_stack(features)
        
        # Standardize if requested
        if standardize:
            scaler = StandardScaler()
            operator_features = scaler.fit_transform(operator_features)
        
        # Store in AnnData
        self.adata.obsm[obsm_key] = operator_features
        
        print(f"Operator features prepared: {operator_features.shape}")
        print(f"Features: {metrics_keys}")
        
        return operator_features
    
    def cluster_operator_space(
        self,
        use_rep: str = 'X_operator',
        n_neighbors: int = 15,
        resolution: float = 1.0,
        key_added: str = 'operator_clusters',
        random_state: int = 0
    ) -> np.ndarray:
        """
        Cluster cells using only operator features.
        
        This demonstrates clustering in operator space, which captures dynamical
        properties rather than expression patterns.
        
        Args:
            use_rep: Key in adata.obsm for operator features
            n_neighbors: Number of neighbors for graph construction
            resolution: Resolution parameter for Leiden clustering
            key_added: Key to store cluster labels in adata.obs
            random_state: Random seed for reproducibility
            
        Returns:
            cluster_labels: Array of cluster labels
        """
        if use_rep not in self.adata.obsm:
            raise ValueError(f"Representation '{use_rep}' not found in adata.obsm. "
                           f"Run prepare_operator_features() first.")
        
        # Build neighborhood graph
        sc.pp.neighbors(
            self.adata,
            use_rep=use_rep,
            n_neighbors=n_neighbors,
            random_state=random_state
        )
        
        # Leiden clustering
        sc.tl.leiden(
            self.adata,
            resolution=resolution,
            key_added=key_added,
            random_state=random_state
        )
        
        n_clusters = self.adata.obs[key_added].nunique()
        
        print(f"Operator-based clustering complete: {key_added}")
        print(f"Number of clusters: {n_clusters}")
        
        return self.adata.obs[key_added].values
    
    def cluster_joint_space(
        self,
        expression_rep: str = 'X_pca',
        operator_rep: str = 'X_operator',
        alpha: float = 0.5,
        n_neighbors: int = 15,
        resolution: float = 1.0,
        key_added: str = 'joint_clusters',
        random_state: int = 0
    ) -> np.ndarray:
        """
        Cluster using combined expression and operator features.
        
        The joint space combines expression (what genes are expressed) with
        operator properties (how stable/plastic the state is).
        
        Args:
            expression_rep: Key in adata.obsm for expression features (e.g., 'X_pca')
            operator_rep: Key in adata.obsm for operator features
            alpha: Weight for operator features (0=expression only, 1=operator only)
            n_neighbors: Number of neighbors
            resolution: Resolution parameter
            key_added: Key to store cluster labels
            random_state: Random seed
            
        Returns:
            cluster_labels: Array of cluster labels
        """
        if expression_rep not in self.adata.obsm:
            raise ValueError(f"Expression representation '{expression_rep}' not found")
        if operator_rep not in self.adata.obsm:
            raise ValueError(f"Operator representation '{operator_rep}' not found")
        
        # Get features
        expr_features = self.adata.obsm[expression_rep]
        oper_features = self.adata.obsm[operator_rep]
        
        # Standardize both
        expr_scaled = StandardScaler().fit_transform(expr_features)
        oper_scaled = StandardScaler().fit_transform(oper_features)
        
        # Weighted combination
        joint_features = (1 - alpha) * expr_scaled + alpha * oper_scaled
        
        # Store
        self.adata.obsm['X_joint'] = joint_features
        
        # Cluster
        sc.pp.neighbors(
            self.adata,
            use_rep='X_joint',
            n_neighbors=n_neighbors,
            random_state=random_state
        )
        sc.tl.leiden(
            self.adata,
            resolution=resolution,
            key_added=key_added,
            random_state=random_state
        )
        
        n_clusters = self.adata.obs[key_added].nunique()
        
        print(f"Joint clustering complete: {key_added}")
        print(f"Alpha (operator weight): {alpha}")
        print(f"Number of clusters: {n_clusters}")
        
        return self.adata.obs[key_added].values
    
    def compare_clustering_quality(
        self,
        methods: Dict[str, str],
        celltype_key: str,
        compute_silhouette: bool = True
    ) -> Dict:
        """
        Compare clustering quality across methods.
        
        This is critical for demonstrating that operator-based clustering
        provides value beyond expression-based clustering.
        
        Args:
            methods: Dict mapping method name to cluster key in adata.obs
                    e.g., {'expression': 'leiden', 'operator': 'operator_clusters'}
            celltype_key: Ground truth cell type labels in adata.obs
            compute_silhouette: Whether to compute silhouette scores (can be slow)
            
        Returns:
            Dictionary with quality metrics for each method:
                - 'ari': Adjusted Rand Index
                - 'nmi': Normalized Mutual Information
                - 'silhouette': Silhouette score (if computed)
                - 'n_clusters': Number of clusters
        """
        if celltype_key not in self.adata.obs:
            raise ValueError(f"Cell type key '{celltype_key}' not found in adata.obs")
        
        true_labels = self.adata.obs[celltype_key].values
        
        results = {}
        print("\n=== Clustering Quality Comparison ===")
        print(f"Ground truth: {celltype_key}")
        print(f"Number of cell types: {len(np.unique(true_labels))}\n")
        
        for method_name, cluster_key in methods.items():
            if cluster_key not in self.adata.obs:
                warnings.warn(f"Cluster key '{cluster_key}' not found, skipping {method_name}")
                continue
            
            pred_labels = self.adata.obs[cluster_key].values
            
            # Compute metrics
            ari = adjusted_rand_score(true_labels, pred_labels)
            nmi = normalized_mutual_info_score(true_labels, pred_labels)
            
            # Silhouette score (requires embedding)
            sil = None
            if compute_silhouette:
                # Determine which embedding to use
                if 'operator' in method_name and 'X_operator' in self.adata.obsm:
                    features = self.adata.obsm['X_operator']
                elif 'joint' in method_name and 'X_joint' in self.adata.obsm:
                    features = self.adata.obsm['X_joint']
                elif 'X_pca' in self.adata.obsm:
                    features = self.adata.obsm['X_pca']
                else:
                    warnings.warn(f"No suitable embedding found for {method_name}, skipping silhouette")
                    features = None
                
                if features is not None:
                    try:
                        sil = silhouette_score(features, pred_labels)
                    except Exception as e:
                        warnings.warn(f"Could not compute silhouette for {method_name}: {e}")
            
            results[method_name] = {
                'ari': ari,
                'nmi': nmi,
                'silhouette': sil,
                'n_clusters': len(np.unique(pred_labels))
            }
            
            # Print results
            print(f"{method_name:15s}:")
            print(f"  ARI:        {ari:.3f}")
            print(f"  NMI:        {nmi:.3f}")
            if sil is not None:
                print(f"  Silhouette: {sil:.3f}")
            print(f"  N clusters: {len(np.unique(pred_labels))}")
            print()
        
        return results
    
    def grid_search_alpha(
        self,
        expression_rep: str = 'X_pca',
        operator_rep: str = 'X_operator',
        alphas: List[float] = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        celltype_key: str = 'cell_type',
        n_neighbors: int = 15,
        resolution: float = 1.0,
        random_state: int = 0
    ) -> Tuple[Dict, float]:
        """
        Grid search over alpha values to find optimal expression/operator balance.
        
        Args:
            expression_rep: Expression representation key
            operator_rep: Operator representation key
            alphas: List of alpha values to test
            celltype_key: Ground truth cell type labels
            n_neighbors: Number of neighbors
            resolution: Clustering resolution
            random_state: Random seed
            
        Returns:
            results: Dictionary with metrics for each alpha
            best_alpha: Alpha value with highest ARI
        """
        results = {}
        best_ari = -1
        best_alpha = 0.0
        
        print("\n=== Grid Search Over Alpha ===")
        print(f"Testing {len(alphas)} alpha values...\n")
        
        for alpha in alphas:
            # Cluster with this alpha
            key = f'joint_alpha_{alpha:.1f}'
            self.cluster_joint_space(
                expression_rep=expression_rep,
                operator_rep=operator_rep,
                alpha=alpha,
                n_neighbors=n_neighbors,
                resolution=resolution,
                key_added=key,
                random_state=random_state
            )
            
            # Evaluate
            true_labels = self.adata.obs[celltype_key].values
            pred_labels = self.adata.obs[key].values
            
            ari = adjusted_rand_score(true_labels, pred_labels)
            nmi = normalized_mutual_info_score(true_labels, pred_labels)
            
            results[alpha] = {
                'ari': ari,
                'nmi': nmi,
                'n_clusters': len(np.unique(pred_labels))
            }
            
            print(f"Alpha={alpha:.1f}: ARI={ari:.3f}, NMI={nmi:.3f}, "
                  f"n_clusters={len(np.unique(pred_labels))}")
            
            if ari > best_ari:
                best_ari = ari
                best_alpha = alpha
        
        print(f"\nBest alpha: {best_alpha:.1f} (ARI={best_ari:.3f})")
        
        return results, best_alpha
    
    def compute_cluster_purity(
        self,
        cluster_key: str,
        celltype_key: str
    ) -> Dict:
        """
        Compute purity of clusters with respect to cell types.
        
        Purity measures how homogeneous each cluster is in terms of cell types.
        
        Args:
            cluster_key: Cluster labels in adata.obs
            celltype_key: Cell type labels in adata.obs
            
        Returns:
            Dictionary with purity metrics for each cluster
        """
        if cluster_key not in self.adata.obs:
            raise ValueError(f"Cluster key '{cluster_key}' not found")
        if celltype_key not in self.adata.obs:
            raise ValueError(f"Cell type key '{celltype_key}' not found")
        
        clusters = self.adata.obs[cluster_key].values
        celltypes = self.adata.obs[celltype_key].values
        
        purity_results = {}
        
        for cluster in np.unique(clusters):
            cluster_mask = clusters == cluster
            cluster_celltypes = celltypes[cluster_mask]
            
            # Most common cell type in this cluster
            unique, counts = np.unique(cluster_celltypes, return_counts=True)
            dominant_celltype = unique[np.argmax(counts)]
            purity = counts.max() / len(cluster_celltypes)
            
            purity_results[cluster] = {
                'size': len(cluster_celltypes),
                'dominant_celltype': dominant_celltype,
                'purity': purity,
                'composition': dict(zip(unique, counts))
            }
        
        # Overall purity
        overall_purity = np.mean([r['purity'] for r in purity_results.values()])
        
        print(f"\n=== Cluster Purity Analysis ===")
        print(f"Overall purity: {overall_purity:.3f}")
        print(f"\nPer-cluster purity:")
        for cluster, stats in purity_results.items():
            print(f"  Cluster {cluster}: {stats['purity']:.3f} "
                  f"(n={stats['size']}, dominant={stats['dominant_celltype']})")
        
        return {
            'per_cluster': purity_results,
            'overall_purity': overall_purity
        }


def quick_operator_clustering(
    adata: ad.AnnData,
    method: str = 'operator',
    alpha: float = 0.5,
    resolution: float = 1.0
) -> ad.AnnData:
    """
    Convenience function for quick operator-based clustering.
    
    Args:
        adata: AnnData object with operator metrics
        method: 'operator', 'expression', or 'joint'
        alpha: Weight for operator features (only for 'joint' method)
        resolution: Clustering resolution
        
    Returns:
        adata: AnnData with clustering results
        
    Example:
        >>> adata = quick_operator_clustering(adata, method='joint', alpha=0.5)
        >>> sc.pl.umap(adata, color='joint_clusters')
    """
    clusterer = OperatorClustering(adata)
    
    if method == 'operator':
        clusterer.prepare_operator_features()
        clusterer.cluster_operator_space(resolution=resolution)
    
    elif method == 'joint':
        clusterer.prepare_operator_features()
        clusterer.cluster_joint_space(alpha=alpha, resolution=resolution)
    
    elif method == 'expression':
        # Standard expression-based clustering
        if 'X_pca' not in adata.obsm:
            raise ValueError("PCA representation not found. Run sc.tl.pca() first.")
        sc.pp.neighbors(adata, use_rep='X_pca')
        sc.tl.leiden(adata, resolution=resolution, key_added='expression_clusters')
    
    else:
        raise ValueError(f"Unknown method: {method}")
    
    return adata
