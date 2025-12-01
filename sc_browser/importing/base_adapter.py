from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Dict, Any

from sc_browser.core.dataset import Dataset

class BaseConfigAdapter(ABC):
    """
    Takes a raw config entry (from JSON) and returns a dict compatible with
    {@link Dataset.from_config_entry()}
    """
    id: str

    @abstractmethod
    def can_handle(self, entry: Dict[str, Any]) -> bool:
        """
        Returns true if this adapter can handle this config entry
        """
        raise NotImplementedError


    @abstractmethod
    def build_dataset(self, entry: Dict[str, Any]) -> Dataset:
        """
        Builds the dataset from the config entry
        """
        raise NotImplementedError



# second iteration, converts .h5ad files

# third iteration - legacy data formats (seurat) converts, etc...


class SimpleSingleCellAdapter(BaseConfigAdapter):
    """
    Very simple adapter for your current single-cell config schema, e.g.:

    {
      "name": "Demo dataset",
      "group": "Example",
      "file": "data/demo.h5ad",
      "cluster_key": "cluster",
      "condition_key": "condition",
      "embedding_key": "X_umap"
    }
    """

    id = "singlecell_simple"

    def can_handle(self, entry: Dict[str, Any]) -> bool:
        required = {"name", "group", "cluster_key", "condition_key", "embedding_key"}
        return required.issubset(entry.keys())

    def build_dataset(self, entry: Dict[str, Any]) -> Dataset:
        # Delegate actual construction to Dataset.from_config_entry
        return Dataset.from_config_entry(entry)