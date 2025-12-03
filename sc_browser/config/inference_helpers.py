from __future__ import annotations

import anndata as AnnData
import pandas as pd


# ---- inference helpers ----

def _pick_obs_column(
    obs: pd.DataFrame,
    candidates: list[str],
    max_unique: int | None = None,
) -> str | None:
    """
    Return the first column in `candidates` that exists in obs, with optional
    limit on cardinality.
    """
    for name in candidates:
        if name in obs.columns:
            if max_unique is not None and obs[name].nunique() > max_unique:
                continue
            return name
    return None


def _infer_cluster_key(adata: AnnData) -> str | None:
    obs = adata.obs
    # common cluster labels from Scanpy/Seurat/pipelines
    candidates = [
        "leiden",
        "louvain",
        "clusters",
        "cluster",
        "seurat_clusters",
        "label"
    ]
    return _pick_obs_column(obs, candidates, max_unique=max(adata.n_obs // 2, 50))


def _infer_condition_key(adata: AnnData) -> str | None:
    obs = adata.obs
    candidates = [
        "condition",
        "cond",
        "treatment",
        "group",
        "status",
    ]
    return _pick_obs_column(obs, candidates, max_unique=50)


def _infer_sample_key(adata: AnnData) -> str | None:
    obs = adata.obs
    candidates = [
        "sample",
        "sample_id",
        "dataset",
        "donor",
        "patient",
    ]
    return _pick_obs_column(obs, candidates)


def _infer_cell_type_key(adata: AnnData) -> str | None:
    obs = adata.obs
    candidates = [
        "cell_type",
        "celltype",
        "celltype_panglao",
        "annotation",
        "cell_identity",
    ]
    return _pick_obs_column(obs, candidates)


def _infer_embedding_key(adata: AnnData) -> str | None:
    """
    Try to pick a sensible .obsm key for embedding.

    Preference:
      - 'X_umap'
      - keys containing 'umap'
      - keys containing 'tsne'
      - otherwise, first key with 2 or 3 dimensions
    """
    obsm_keys = list(adata.obsm.keys())
    if not obsm_keys:
        return None

    if "X_umap" in obsm_keys:
        return "X_umap"

    # case-insensitive search
    for key in obsm_keys:
        if "umap" in key.lower():
            return key

    for key in obsm_keys:
        if "tsne" in key.lower():
            return key

    # fallback: any 2D/3D embedding
    for key in obsm_keys:
        arr = adata.obsm[key]
        if hasattr(arr, "shape") and arr.shape[1] in (2, 3):
            return key

    return obsm_keys[0]