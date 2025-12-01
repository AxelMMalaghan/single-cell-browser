from __future__ import annotations

from typing import Dict, Any

from sc_browser.core.dataset import Dataset
from sc_browser.core.base_adapter import BaseConfigAdapter


class SimpleSingleCellAdapter(BaseConfigAdapter):
    """
    Very simple adapter for the current single-cell JSON config schema.

    Expected per-dataset entry (from config.json):

        {
          "name": "Demo dataset",
          "group": "Example",
          "file": "data/demo.h5ad",        # or "file_path"
          "cluster_key": "cluster",
          "condition_key": "condition",
          "embedding_key": "X_umap"
        }

    Responsibilities:
    - Decide if it can handle this entry (can_handle)
    - Map that entry into a core Dataset (build_dataset)
    """

    id = "singlecell_simple"

    def can_handle(self, entry: Dict[str, Any]) -> bool:
        """
        Return True if this entry looks like our v1 single-cell schema.
        We don't do strict validation here, just enough to pick the adapter.
        """
        return entry.get("schema") == "anndata_mapped"

    def build_dataset(self, entry: Dict[str, Any]) -> Dataset:
        """
        Build a Dataset from this config entry.

        All the actual I/O and AnnData wiring is delegated to Dataset.from_config_entry,
        so this adapter only decides *that* it can handle the schema.
        """
        return Dataset.from_config_entry(entry)