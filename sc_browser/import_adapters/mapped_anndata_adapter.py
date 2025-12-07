from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

from sc_browser.core.base_adapter import BaseConfigAdapter
from sc_browser.config.model import DatasetConfig
from sc_browser.core.dataset import Dataset
from sc_browser.core.dataset_loader import from_config


class MappedAnnDataAdapter(BaseConfigAdapter):
    """
    Adapter for config entries that explicitly describe an AnnData dataset
    using the 'anndata_mapped' schema.

    This is essentially a thin wrapper around DatasetConfig + from_config.
    """

    id = "anndata_mapped"

    def can_handle(self, entry: Dict[str, Any]) -> bool:
        """
        Handle entries of the form:
        {
          "schema": "anndata_mapped",
          "name": "...",
          "file" | "path" | "file_path": "path/to/file.h5ad",
          "obs_columns": { ... },
          "group": "...",
          "embedding_key": "X_umap"
        }
        """
        return (
            entry.get("schema") == self.id
            and ("file" in entry or "path" in entry or "file_path" in entry)
        )

    def build_dataset(self, entry: Dict[str, Any]) -> Dataset:
        """
        Build a Dataset from a raw config dict using DatasetConfig + from_config.
        """
        source_path = Path(
            entry.get("file")
            or entry.get("path")
            or entry.get("file_path", ".")
        )
        cfg = DatasetConfig(
            raw=entry,
            source_path=source_path,
            index=0,  # if used in batch configs, caller can inject proper index
        )
        return from_config(cfg)
