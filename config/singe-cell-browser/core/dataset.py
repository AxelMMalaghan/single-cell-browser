from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

import anndata as ad
import numpy as np
import pandas as pd

@dataclass
class Dataset:
    """
    Wraps an Anndata Object plus configuration metadata.
    Represents one experiment
    """
    name: str
    group: str
    adata: ad.AnnData
    config: Dict[str, Any]

    @classmethod
    def from_config_entry(cls, entry: Dict[str, Any]) -> "Dataset":
        """
        Create from a JSON config block
        """
        file_path = Path(entry["file_path"])
        adata = ad.read_h5ad(file_path)
        return cls(
            name=entry["name"],
            group=entry.get("group", "default"),
            adata=adata,
            config=entry
        )


    @property
    def genes(self) -> List[str]:
       return list(self.adata.var_names)


    @property
    def cluster_key(self) -> str:
        return self.config.get(["cluster_key", "cluster"])

    @property
    def condition_key(self) -> str:
        return self.config.get(["condition_key", "condition"])

    @property
    def embedding_key(self) -> str:
        return self.config.get(["embedding_key", "embedding"])

    @property
    def clusters(self) -> pd.Series:
        return self.adata.obs[self.cluster_key].astype(str)

    @property
    def conditions(self) -> pd.Series:
        return self.adata.obs[self.condition_key].astype(str)

    def subset(self, clusters:Optional[List[str]] = None, conditions:Optional[List[str]]=None) -> "Dataset":
        """
        Return a view of the dataset filtered by cluster/condition
        All views should go through this, not re-implment filtering
        """

        mask = np.ones(self.adata.n_obs, dtype=bool)

        if clusters:
            mask &= self.clusters.isin(clusters).to_numpy()

        if conditions:
            mask &= self.conditions.isin(conditions).to_numpy()

        sub_adata = self.adata[mask].copy()

        return Dataset(name=self.name, group=self.group, adata=sub_adata, config=self.config)

