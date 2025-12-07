from __future__ import annotations

from typing import Optional, List, Dict

import anndata
import pandas as pd

from sc_browser.config.model import DatasetConfig, ObsColumns

class AnnDataAdapter:
    """
    Wraps access to a .h5ad AnnData file based on semantic column mappings.
    """

    def __init__(self, config: DatasetConfig):
        self.config = config
        self.adata = anndata.read_h5ad(config.path)
        self.obs_columns: ObsColumns = config.obs_columns

    def get_obs_column(self, role: str) -> Optional[pd.Series]:
        col_name = getattr(self.obs_columns, role, None)
        if col_name and col_name in self.adata.obs:
            return self.adata.obs[col_name]
        return None

    def get_cell_ids(self) -> pd.Index:
        return self.adata.obs_names

    def get_cluster_labels(self) -> Optional[pd.Series]:
        return self.get_obs_column("cluster")

    def get_condition_labels(self) -> Optional[pd.Series]:
        return self.get_obs_column("condition")

    def get_sample_ids(self) -> Optional[pd.Series]:
        return self.get_obs_column("sample")

    def get_cell_types(self) -> Optional[pd.Series]:
        return self.get_obs_column("cell_type")

    def get_embedding_keys(self) -> List[str]:
        return list(self.adata.obsm.keys())

    def get_embedding(self, key: str) -> Optional[pd.DataFrame]:
        if key in self.adata.obsm:
            emb = self.adata.obsm[key]
            return pd.DataFrame(emb, index=self.adata.obs_names)
        return None

    def summary(self) -> Dict[str, any]:
        return {
            "name": self.config.name,
            "n_obs": self.adata.n_obs,
            "n_vars": self.adata.n_vars,
            "obs_columns": list(self.adata.obs.columns),
            "obsm_keys": list(self.adata.obsm.keys()),
        }
