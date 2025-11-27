# sc_browser/core/dataset.py

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import anndata as ad
import pandas as pd


class Dataset:
    """
    Domain wrapper around an AnnData object, plus config metadata.

    This class hides all the AnnData- and schema specific details from the rest of the app

    Design goals:
    - Provide a stale, typed interface for views and other users to access clusters, conditions, genes, and embeddings
    without depending on the underlying AnnData structure.
    - Centralise config-driven concerns such as column keys so that changing schema only requires changes in one place.
    - Cheap creation of filtered views of the data via subset() without exposing raw AnnData objects.
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

    @classmethod
    def from_config_entry(cls, entry: dict) -> "Dataset":
        """
        Build a Dataset from a config.json file.

        The entry is expected to have at least:
        - "name":           display name for dataset
        - "group":          grouping label for dataset
        - "cluster_key":    obs column name for clusters
        - "condition_key":  obs column name for conditions
        - "embedding_key":   obsm key for embedding

        It also needs a file path to an .h5ad file, given by either:
        - "file_path"
        - "file"

        This factory keeps all file I/O and schema wiring in one place so that callers (and views) just deal with the Dataset Abstraction

        Raises:
            KeyError: if file or file path does not exist
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
            embedding_key=entry["embedding_key"]
        )


    def subset(
        self,
        clusters: Optional[List[str]] = None,
        conditions: Optional[List[str]] = None,
    ) -> "Dataset":
        """
        Return a new Dataset filtered by cluster and/or condition.

        This method builds a boolean mask over obs using the configured cluster and condition keys, subsets the
        underlying AnnData and wraps the result in a new Dataset instances the preserves all metadata.
        (name, group, key names, embedding key)

        :param clusters: cluster labels as a string or list of cluster labels as a string, if None or empty, no filtering is applied
        :param conditions: condition labels as a string or list of condition labels as a string, if None or empty, no filtering is applied

        :return a new Dataset object backed by a copy of the filtered AnnData

        Notes:
        This keeps callers working with the Dataset abstraction rather than slicing AnnData directly
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


    # --- Getters ---

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
            columns=[f"dim{i + 1}" for i in range(emb.shape[1])],
        )