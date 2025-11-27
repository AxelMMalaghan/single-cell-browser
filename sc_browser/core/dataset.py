# sc_browser/core/dataset.py

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import anndata as ad
import pandas as pd


class Dataset:
    """
    Domain wrapper around an AnnData object, plus config metadata.

    This isolates all the "what is the cluster column called", "which embedding
    do we use", etc. so the rest of the app doesn't depend directly on AnnData
    internals.
    """

    def __init__(
        self,
        name: str,
        group: str,
        adata: ad.AnnData,
        cluster_key: str,
        condition_key: str,
        embedding_key: str,
    ) -> None:
        self.name = name
        self.group = group
        self.adata = adata
        self.cluster_key = cluster_key
        self.condition_key = condition_key
        self.embedding_key = embedding_key

    # ---------- factory from config ----------

    @classmethod
    def from_config_entry(cls, entry: dict) -> "Dataset":
        """
        Build a Dataset from a config.json entry.

        Supports either "file" or "file_path" as the path key.
        """
        file_value = entry.get("file_path") or entry.get("file")
        if file_value is None:
            raise KeyError("Config entry must contain 'file' or 'file_path'")

        file_path = Path(file_value)
        adata = ad.read_h5ad(file_path)

        return cls(
            name=entry["name"],
            group=entry["group"],
            adata=adata,
            cluster_key=entry["cluster_key"],
            condition_key=entry["condition_key"],
            embedding_key=entry["embedding_key"],
        )

    # ---------- convenience accessors ----------

    @property
    def clusters(self) -> pd.Series:
        """Return cluster labels as a string Series."""
        return self.adata.obs[self.cluster_key].astype(str)

    @property
    def conditions(self) -> pd.Series:
        """Return condition labels as a string Series."""
        return self.adata.obs[self.condition_key].astype(str)

    @property
    def genes(self) -> pd.Index:
        """Return gene names for this a given dataset"""
        return self.adata.var_names

    @property
    def embedding(self) -> pd.DataFrame:
        """Return the 2D embedding as a DataFrame with dim1/dim2 columns."""
        emb = self.adata.obsm[self.embedding_key]
        return pd.DataFrame(
            emb,
            index=self.adata.obs.index,
            columns=[f"dim{i+1}" for i in range(emb.shape[1])],
        )

    # ---------- subsetting ----------

    def subset(
        self,
        clusters: Optional[List[str]] = None,
        conditions: Optional[List[str]] = None,
    ) -> "Dataset":
        """
        Return a new Dataset filtered by cluster and/or condition.
        """
        mask = pd.Series(True, index=self.adata.obs.index)

        if clusters:
            mask = mask & self.clusters.isin(clusters)
        if conditions:
            mask = mask & self.conditions.isin(conditions)

        sub = self.adata[mask].copy()

        return Dataset(
            name=self.name,
            group=self.group,
            adata=sub,
            cluster_key=self.cluster_key,
            condition_key=self.condition_key,
            embedding_key=self.embedding_key,
        )
