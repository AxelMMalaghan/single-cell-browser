from __future__ import annotations
from pathlib import Path
from typing import Dict, Optional, List

import anndata as ad
import pandas as pd


class Dataset:

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


    def subset(
            self,
            clusters: Optional[List[str]] = None,
            conditions: Optional[List[str]] = None,
            samples: Optional[List[str]] = None,
            cell_types: Optional[List[str]] = None,
    ) -> "Dataset":
        """
        Return a new Dataset filtered by cluster, condition, sample, and/or cell type.

        This method builds a boolean mask over `.obs` using the configured keys,
        subsets the underlying AnnData, and wraps the result in a new Dataset instance.
        """
        obs = self.adata.obs
        mask = pd.Series(True, index=obs.index)

        # Base filters using standard keys
        if clusters:
            mask &= obs[self.cluster_key].astype(str).isin(clusters)

        if conditions:
            mask &= obs[self.condition_key].astype(str).isin(conditions)


        # Optional filters based on configured obs_columns (dict of semantic name -> obs column)
        sample_key = self.obs_columns.get("sample")
        if samples and sample_key and sample_key in obs.columns:
            mask &= obs[sample_key].astype(str).isin(samples)

        cell_type_key = self.obs_columns.get("cell_type")
        if cell_types and cell_type_key and cell_type_key in obs.columns:
            mask &= obs[cell_type_key].astype(str).isin(cell_types)

        # Subset the data
        sub = self.adata[mask].copy()

        return Dataset(
            name=self.name,
            group=self.group,
            adata=sub,
            cluster_key=self.cluster_key,
            condition_key=self.condition_key,
            embedding_key=self.embedding_key,
            obs_columns=self.obs_columns,  # preserve extended keys
            file_path=self.file_path,
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

    def get_obs_column(self, key: str) -> Optional[str]:
        return self.obs_columns.get(key)

