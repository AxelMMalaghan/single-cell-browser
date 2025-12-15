from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, TYPE_CHECKING

import anndata as ad
import numpy as np
import pandas as pd

from sc_browser.config.model import ObsColumns

if TYPE_CHECKING:
    from sc_browser.core.filter_state import FilterState


@dataclass(frozen=True)
class ValidSets:
    """
    Cached valid values for UI validation / sanitisation.
    Computing .unique() on large Series can be expensive, so we do it once per Dataset.
    """
    clusters: set[str]
    conditions: set[str]
    samples: set[str]
    cell_types: set[str]
    genes: set[str]
    embeddings: set[str]


class Dataset:
    """
    Unified dataset abstraction used throughout the browser.

    Includes:
    - Standardised access to .obs columns
    - Efficient cached subsetting of AnnData
    - Cached expression extraction (cells × genes)
    - Sparse-aware handling where possible
    """

    MAX_SUBSET_CACHE = 128
    MAX_EXPR_CACHE = 128

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
        obs_columns: Optional[Any] = None,
        file_path: Optional[Path] = None,
    ) -> None:
        self.name = name
        self.group = group
        self.adata = adata
        self.cluster_key = cluster_key
        self.condition_key = condition_key
        self.embedding_key = embedding_key
        self.file_path = file_path

        # ---------------------------------------------------------------------
        # Normalise obs_columns to a dict[str, str]
        # Accept both:
        #   - dict-like: {"sample": "sample_col", "cell_type": "cell_type_col", ...}
        #   - ObsColumns dataclass instance
        # ---------------------------------------------------------------------
        self.obs_columns: Dict[str, str] = self._normalise_obs_columns(obs_columns)

        # ---------------------------------------------------------------------
        # Pre-normalise all relevant .obs columns to avoid repeated astype(str)
        # ---------------------------------------------------------------------
        obs = adata.obs

        self._cluster_series: Optional[pd.Series] = (
            obs[cluster_key].astype(str).copy()
            if cluster_key is not None and cluster_key in obs.columns
            else None
        )

        self._condition_series: Optional[pd.Series] = (
            obs[condition_key].astype(str).copy()
            if condition_key is not None and condition_key in obs.columns
            else None
        )

        sample_key = self.obs_columns.get("sample")
        self._sample_series: Optional[pd.Series] = (
            obs[sample_key].astype(str).copy()
            if sample_key and sample_key in obs.columns
            else None
        )

        cell_type_key = self.obs_columns.get("cell_type")
        self._cell_type_series: Optional[pd.Series] = (
            obs[cell_type_key].astype(str).copy()
            if cell_type_key and cell_type_key in obs.columns
            else None
        )

        # ---------------------------------------------------------------------
        # Caches
        # ---------------------------------------------------------------------
        # Cache of filtered Dataset objects
        self._subset_cache: Dict[
            Tuple[Tuple[str, ...], Tuple[str, ...], Tuple[str, ...], Tuple[str, ...]],
            "Dataset",
        ] = {}

        # Cache for expression matrices (per subset, per gene set)
        self._expr_cache: Dict[Tuple[str, ...], pd.DataFrame] = {}

        # Cached valid values for UI sanitisation
        self._valid_sets: Optional[ValidSets] = None

    # -------------------------------------------------------------------------
    # Internal: normalise obs_columns
    # -------------------------------------------------------------------------
    def _normalise_obs_columns(self, obs_columns: Optional[Any]) -> Dict[str, str]:
        """
        Accept either:
        - dict-like: {"sample": "...", "cell_type": "...", ...}
        - ObsColumns dataclass instance

        and always return a flat dict[str, str] with only non-None entries.
        """
        if obs_columns is None:
            return {}

        if isinstance(obs_columns, dict):
            return {k: v for k, v in obs_columns.items() if v is not None}

        if isinstance(obs_columns, ObsColumns):
            mapping: Dict[str, str] = {}
            if obs_columns.sample is not None:
                mapping["sample"] = obs_columns.sample
            if obs_columns.cell_type is not None:
                mapping["cell_type"] = obs_columns.cell_type
            if obs_columns.batch is not None:
                mapping["batch"] = obs_columns.batch
            if obs_columns.cell_id is not None:
                mapping["cell_id"] = obs_columns.cell_id
            if obs_columns.cluster is not None:
                mapping["cluster"] = obs_columns.cluster
            if obs_columns.condition is not None:
                mapping["condition"] = obs_columns.condition
            return mapping

        try:
            return {k: v for k, v in dict(obs_columns).items() if v is not None}
        except Exception:
            return {}

    # -------------------------------------------------------------------------
    # Cached valid values for UI sanitisation
    # -------------------------------------------------------------------------
    def valid_sets(self) -> ValidSets:
        """
        Return cached valid values for dropdown sanitisation.

        This avoids recomputing `.unique()` in callbacks on every UI change.
        """
        if self._valid_sets is not None:
            return self._valid_sets

        clusters = set(self._cluster_series.unique()) if self._cluster_series is not None else set()
        conditions = set(self._condition_series.unique()) if self._condition_series is not None else set()
        samples = set(self._sample_series.unique()) if self._sample_series is not None else set()
        cell_types = set(self._cell_type_series.unique()) if self._cell_type_series is not None else set()
        genes = set(map(str, self.adata.var_names))
        embeddings = set(map(str, self.adata.obsm.keys()))

        self._valid_sets = ValidSets(
            clusters=clusters,
            conditions=conditions,
            samples=samples,
            cell_types=cell_types,
            genes=genes,
            embeddings=embeddings,
        )
        return self._valid_sets

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
        *,
        copy_adata: bool = False,
    ) -> "Dataset":
        """
        Efficiently return a new Dataset filtered by cluster, condition, sample,
        and/or cell type. Subsets are cached to avoid recomputation.

        Args:
            copy_adata: If True, forces `adata[mask].copy()` (safer if any code mutates AnnData).
                       Default False keeps views (lower memory) assuming read-only usage.
        """
        key = self._subset_cache_key(clusters, conditions, samples, cell_types)
        cached = self._subset_cache.get(key)
        if cached is not None:
            return cached

        adata = self.adata
        n = adata.n_obs

        mask = np.ones(n, dtype=bool)

        if clusters and self._cluster_series is not None:
            mask &= self._cluster_series.isin(clusters).to_numpy()

        if conditions and self._condition_series is not None:
            mask &= self._condition_series.isin(conditions).to_numpy()

        if samples and self._sample_series is not None:
            mask &= self._sample_series.isin(samples).to_numpy()

        if cell_types and self._cell_type_series is not None:
            mask &= self._cell_type_series.isin(cell_types).to_numpy()

        sub_adata = adata[mask].copy() if copy_adata else adata[mask]

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

        self._subset_cache[key] = subset_ds

        # Prevent unbounded growth
        if len(self._subset_cache) > self.MAX_SUBSET_CACHE:
            self._subset_cache.clear()

        return subset_ds

    # -------------------------------------------------------------------------
    # Centralised subsetting based on FilterState
    # -------------------------------------------------------------------------
    def subset_for_state(self, state: "FilterState") -> "Dataset":
        """
        Convenience wrapper to subset this Dataset based on a FilterState.
        """
        cell_types = getattr(state, "cell_types", None) or getattr(state, "celltypes", None)

        return self.subset(
            clusters=getattr(state, "clusters", None) or None,
            conditions=getattr(state, "conditions", None) or None,
            samples=getattr(state, "samples", None) or None,
            cell_types=cell_types or None,
        )

    # -------------------------------------------------------------------------
    # Expression Matrix Extraction (cached)
    # -------------------------------------------------------------------------
    def expression_matrix(self, genes: Sequence[str]) -> pd.DataFrame:
        """
        Return an expression matrix (cells × genes) for the given genes.

        Notes:
        - Cache key is order-insensitive (genes are normalised to a sorted unique tuple).
        - If none of the genes exist, returns an empty DataFrame indexed by obs_names.
        """
        key = tuple(sorted(set(map(str, genes))))
        if key in self._expr_cache:
            return self._expr_cache[key]

        adata = self.adata

        var_mask = adata.var_names.isin(list(key))
        if not var_mask.any():
            df = pd.DataFrame(index=adata.obs_names)
            self._expr_cache[key] = df
            return df

        selected_genes = adata.var_names[var_mask]
        X = adata[:, var_mask].X  # sparse or dense slice

        if hasattr(X, "toarray"):
            X = X.toarray()

        df = pd.DataFrame(X, index=adata.obs_names, columns=selected_genes)

        self._expr_cache[key] = df

        if len(self._expr_cache) > self.MAX_EXPR_CACHE:
            self._expr_cache.clear()

        return df

    # -------------------------------------------------------------------------
    # Cache management
    # -------------------------------------------------------------------------
    def clear_caches(self) -> None:
        """Reset caches for subsets, expression matrices, and valid value sets."""
        self._subset_cache.clear()
        self._expr_cache.clear()
        self._valid_sets = None

    def refresh_obs_indexes(self) -> None:
        """
        Rebuild cached, string-normalised obs Series after mapping changes.
        Also clears caches and valid_sets.
        """
        obs = self.adata.obs

        self._cluster_series = (
            obs[self.cluster_key].astype(str).copy()
            if self.cluster_key and self.cluster_key in obs.columns
            else None
        )

        self._condition_series = (
            obs[self.condition_key].astype(str).copy()
            if self.condition_key and self.condition_key in obs.columns
            else None
        )

        sample_key = self.obs_columns.get("sample")
        self._sample_series = (
            obs[sample_key].astype(str).copy()
            if sample_key and sample_key in obs.columns
            else None
        )

        cell_type_key = self.obs_columns.get("cell_type")
        self._cell_type_series = (
            obs[cell_type_key].astype(str).copy()
            if cell_type_key and cell_type_key in obs.columns
            else None
        )

        self.clear_caches()

    # -------------------------------------------------------------------------
    # Properties & getters
    # -------------------------------------------------------------------------
    @property
    def clusters(self) -> Optional[pd.Series]:
        """Return cluster labels (string-normalised), if configured."""
        return self._cluster_series

    @property
    def conditions(self) -> Optional[pd.Series]:
        """Return condition labels (string-normalised), if configured."""
        return self._condition_series

    @property
    def samples(self) -> Optional[pd.Series]:
        """Return sample labels (string-normalised), if configured."""
        return self._sample_series

    @property
    def cell_types(self) -> Optional[pd.Series]:
        """Return cell type labels (string-normalised), if configured."""
        return self._cell_type_series

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

        if isinstance(emb, pd.DataFrame):
            df = emb.copy()
            if df.shape[1] > 2:
                df = df.iloc[:, :2]
            df.columns = [f"dim{i + 1}" for i in range(df.shape[1])]
            return df

        arr = np.asarray(emb)
        if arr.ndim != 2:
            raise ValueError(f"Embedding '{emb_key}' must be 2D, got shape {arr.shape}")

        cols = [f"dim{i + 1}" for i in range(arr.shape[1])]
        return pd.DataFrame(arr, index=self.adata.obs.index, columns=cols)

    def get_embedding_matrix(self, key: str) -> np.ndarray:
        """Return embedding matrix as a numpy array."""
        emb = self.adata.obsm[key]
        if isinstance(emb, pd.DataFrame):
            return emb.to_numpy()
        return np.asarray(emb)

    def get_embedding_labels(self, key: str) -> List[str]:
        """Return names of embedding axes."""
        emb = self.adata.obsm[key]

        if isinstance(emb, pd.DataFrame):
            return [str(c) for c in emb.columns]

        arr = np.asarray(emb)
        n = arr.shape[1]
        return [f"Dim {i + 1}" for i in range(n)]

    @property
    def embedding(self) -> pd.DataFrame:
        """Default embedding using this dataset's embedding_key."""
        return self.get_embedding()

    def get_obs_column(self, key: str) -> Optional[str]:
        """
        Backwards-compatible dict-style access to semantic obs columns.
        """
        return self.obs_columns.get(key)
