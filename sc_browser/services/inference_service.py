from __future__ import annotations

from typing import List, Optional

import anndata as ad


def _find_match(columns: List[str], candidates: List[str]) -> Optional[str]:
    """
    Helper to find the first matching column from a list of candidates.
    1. Exact match
    2. Case-insensitive match
    """
    # 1. Exact match
    for cand in candidates:
        if cand in columns:
            return cand

    # 2. Case-insensitive match
    col_map = {c.lower(): c for c in columns}
    for cand in candidates:
        if cand.lower() in col_map:
            return col_map[cand.lower()]

    return None


def infer_cluster_key(adata: ad.AnnData) -> Optional[str]:
    candidates = [
        "leiden",
        "louvain",
        "cluster",
        "clusters",
        "seurat_clusters",
        "RNA_snn_res.0.8",
        "res.0.5",
    ]
    return _find_match(list(adata.obs.columns), candidates)


def infer_condition_key(adata: ad.AnnData) -> Optional[str]:
    candidates = [
        "condition",
        "treatment",
        "group",
        "disease",
        "genotype",
        "status",
        "diagnosis",
        "stim",
    ]
    return _find_match(list(adata.obs.columns), candidates)


def infer_sample_key(adata: ad.AnnData) -> Optional[str]:
    candidates = [
        "sample",
        "sample_id",
        "donor",
        "patient",
        "orig.ident",
        "batch",
        "replicate",
    ]
    return _find_match(list(adata.obs.columns), candidates)


def infer_cell_type_key(adata: ad.AnnData) -> Optional[str]:
    candidates = [
        "cell_type",
        "celltype",
        "annotation",
        "cell_annotation",
        "broad_cell_type",
        "predicted.id",
    ]
    return _find_match(list(adata.obs.columns), candidates)


def infer_embedding_key(adata: ad.AnnData) -> Optional[str]:
    """
    Prefer UMAP > t-SNE > PCA.
    """
    keys = list(adata.obsm.keys())

    # Priority list
    priorities = ["X_umap", "umap", "X_tsne", "tsne", "X_pca", "pca"]

    for p in priorities:
        if p in keys:
            return p

    # Fallback: just take the first one if any exist
    return keys[0] if keys else None
