from __future__ import annotations
from pathlib import Path
from typing import Dict, Optional, List, Tuple, Sequence

import anndata as ad
import numpy as np
import pandas as pd


class Dataset:
    """
    Unified dataset abstraction used throughout the browser.

    Includes:
    - Standardised access to .obs columns
    - Efficient cached subsetting of AnnData
    - Cached expression extraction (cells × genes)
    - Sparse-aware handling where possible
    """

    # -------------------------------------------------------------------------
    # Constructor
    # -------------------------------------------------------------------------
    def __init__(
        self,
        name: str,
        group: str,
        adata: ad.AnnData,
        cluster_key: Optional[str],
        condition_key: Optional[str],
        embedding_key: Optional[str],
        obs_columns: Optional[Dict[str, str]] = None,
        file_path: Optional[Path] = None,
    ) -> None:

        self.name = name
        self.group = group
        self.adata = adata
        self.cluster_key = cluster_key
        self.condition_key = condition_key
        self.embedding_key = embedding_key
        self.obs_columns = obs_columns or {}
        self.file_path = file_path

        # ---------------------------------------------------------------------
        # Pre-normalise all relevant .obs columns to avoid repeated astype(str)
        # ---------------------------------------------------------------------
        obs = adata.obs

        self._cluster_series = (
            obs[cluster_key].astype(str).copy()
            if cluster_key is not None else None
        )

        self._condition_series = (
            obs[condition_key].astype(str).copy()
            if condition_key is not None else None
        )

        sample_key = self.obs_columns.get("sample")
        self._sample_series = (
            obs[sample_key].astype(str).copy()
            if sample_key and sample_key in obs.columns else None
        )

        cell_type_key = self.obs_columns.get("cell_type")
        self._cell_type_series = (
            obs[cell_type_key].astype(str).copy()
            if cell_type_key and cell_type_key in obs.columns else None
        )

        # ---------------------------------------------------------------------
        # Caches
        # ---------------------------------------------------------------------
        # Cache of filtered Dataset objects
        self._subset_cache: Dict[
            Tuple[Tuple[str, ...], Tuple[str, ...], Tuple[str, ...], Tuple[str, ...]],
            "Dataset"
        ] = {}

        # Cache for expression matrices (per subset, per gene set)
        self._expr_cache: Dict[Tuple[str, ...], pd.DataFrame] = {}

    # -------------------------------------------------------------------------
    # Helper: build cache key for subsetting
    # -------------------------------------------------------------------------
    def _subset_cache_key(
        self,
        clusters: Optional[List[str]],
        conditions: Optional[List[str]],
        samples: Optional[List[str]],
        cell_types: Optional[List[str]],
    ) -> Tuple[Tuple[str, ...], Tuple[str, ...], Tuple[str, ...], Tuple[str, ...]]:

        def norm(values: Optional[List[str]]) -> Tuple[str, ...]:
            if not values:
                return ()
            return tuple(sorted(str(v) for v in values))

        return (
            norm(clusters),
            norm(conditions),
            norm(samples),
            norm(cell_types),
        )

    # -------------------------------------------------------------------------
    # Subsetting with caching
    # -------------------------------------------------------------------------
    def subset(
        self,
        clusters: Optional[List[str]] = None,
        conditions: Optional[List[str]] = None,
        samples: Optional[List[str]] = None,
        cell_types: Optional[List[str]] = None,
    ) -> "Dataset":
        """
        Efficiently return a new Dataset filtered by cluster, condition, sample,
        and/or cell type. Subsets are cached to avoid recomputation.
        """
        # Fast-path: return cached subset
        key = self._subset_cache_key(clusters, conditions, samples, cell_types)
        cached = self._subset_cache.get(key)
        if cached is not None:
            return cached

        adata = self.adata
        n = adata.n_obs

        # Start with full mask
        mask = np.ones(n, dtype=bool)

        # Apply pre-normalised filters
        if clusters and self._cluster_series is not None:
            mask &= self._cluster_series.isin(clusters).to_numpy()

        if conditions and self._condition_series is not None:
            mask &= self._condition_series.isin(conditions).to_numpy()

        if samples and self._sample_series is not None:
            mask &= self._sample_series.isin(samples).to_numpy()

        if cell_types and self._cell_type_series is not None:
            mask &= self._cell_type_series.isin(cell_types).to_numpy()

        # Slice WITHOUT copying (views are fine since we never mutate)
        sub_adata = adata[mask]

        subset_ds = Dataset(
            name=self.name,
            group=self.group,
            adata=sub_adata,
            cluster_key=self.cluster_key,
            condition_key=self.condition_key,
            embedding_key=self.embedding_key,
            obs_columns=self.obs_columns,
            file_path=self.file_path,
        )

        # Register cache entry
        self._subset_cache[key] = subset_ds
        return subset_ds

    # -------------------------------------------------------------------------
    # Expression Matrix Extraction (cached)
    # -------------------------------------------------------------------------
    def expression_matrix(self, genes: Sequence[str]) -> pd.DataFrame:
        """
        Return an expression matrix (cells × genes) for the given genes,
        using caching and restricting computation to the relevant subset.
        """
        key = tuple(genes)
        if key in self._expr_cache:
            return self._expr_cache[key]

        adata = self.adata

        # Restrict to genes that actually exist
        var_mask = adata.var_names.isin(genes)
        if not var_mask.any():
            df = pd.DataFrame(index=adata.obs_names)
            self._expr_cache[key] = df
            return df

        selected_genes = adata.var_names[var_mask]
        X = adata[:, var_mask].X  # sparse or dense slice

        # Convert to dense only if needed
        if hasattr(X, "toarray"):
            X = X.toarray()

        df = pd.DataFrame(X, index=adata.obs_names, columns=selected_genes)

        self._expr_cache[key] = df
        return df

    # -------------------------------------------------------------------------
    # Cache management
    # -------------------------------------------------------------------------
    def clear_caches(self) -> None:
        """Reset caches for subsets and expression matrices."""
        self._subset_cache.clear()
        self._expr_cache.clear()

    # -------------------------------------------------------------------------
    # Properties & getters
    # -------------------------------------------------------------------------
    @property
    def clusters(self) -> pd.Series:
        """Return cluster labels."""
        return self._cluster_series

    @property
    def conditions(self) -> pd.Series:
        """Return condition labels."""
        return self._condition_series

    @property
    def genes(self) -> pd.Index:
        """Return gene names."""
        return self.adata.var_names

    def get_embedding(self, key: Optional[str] = None) -> pd.DataFrame:
        """
        Return the embedding for the given obsm key as a DataFrame with dim1/dim2,...

        If key is None, uses this dataset's configured embedding_key.
        """
        emb_key = key or self.embedding_key
        if emb_key is None:
            raise ValueError("No embedding key specified for this dataset")

        if emb_key not in self.adata.obsm:
            raise ValueError(f"Embedding '{emb_key}' not found in adata.obsm")

        emb = self.adata.obsm[emb_key]

        # If already a DataFrame
        if isinstance(emb, pd.DataFrame):
            df = emb.copy()
            if df.shape[1] > 2:
                df = df.iloc[:, :2]
            df.columns = [f"dim{i + 1}" for i in range(df.shape[1])]
            return df

        # Otherwise assume array-like
        arr = np.asarray(emb)
        if arr.ndim != 2:
            raise ValueError(f"Embedding '{emb_key}' must be 2D, got shape {arr.shape}")

        n_dims = arr.shape[1]
        cols = [f"dim{i + 1}" for i in range(n_dims)]
        return pd.DataFrame(arr, index=self.adata.obs.index, columns=cols)

    @property
    def embedding(self) -> pd.DataFrame:
        """Default embedding using this dataset's embedding_key."""
        return self.get_embedding()

    def get_obs_column(self, key: str) -> Optional[str]:
        return self.obs_columns.get(key)
