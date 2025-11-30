# sc_browser/importing/simple_singlecell_adapter.py

from __future__ import annotations

from typing import Any, Mapping, Tuple, List

from ..config_model import GlobalConfig
from ..core.dataset import Dataset
from .base_adapter import BaseConfigAdapter  # whatever you called it


class SimpleSingleCellJsonAdapter(BaseConfigAdapter):
    """
    Adapter for the current 'v1' single-cell JSON config.

    Expected shape of the raw JSON:

    {
      "global": {
        "ui_title": "Single-cell Browser (Mock)",
        "default_group": "Example"
      },
      "datasets": [
        {
          "name": "Demo dataset",
          "group": "Example",
          "file": "data/demo.h5ad",
          "cluster_key": "cluster",
          "condition_key": "condition",
          "embedding_key": "X_umap"
        },
        ...
      ]
    }

    Responsibilities:
    - Decide if it can handle a given raw config (can_handle)
    - Map that config into:
        * GlobalConfig    (UI metadata)
        * list[Dataset]   (domain objects wrapping AnnData)
    """

    adapter_id = "simple_singlecell_json"

    # --------- "Can I handle this config?" ---------

    def can_handle(self, raw: Mapping[str, Any]) -> bool:
        """
        Very lightweight schema check.
        We don't do strict validation here, just enough to select the right adapter.
        """
        if not isinstance(raw, Mapping):
            return False

        datasets = raw.get("datasets")
        if not isinstance(datasets, list) or not datasets:
            return False

        first = datasets[0]
        if not isinstance(first, Mapping):
            return False

        required_keys = {
            "name",
            "group",
            # we accept either file or file_path
            "cluster_key",
            "condition_key",
            "embedding_key",
        }

        # Must at least have the logical keys – file/file_path is checked later
        if not required_keys.issubset(first.keys()):
            return False

        if not ("file" in first or "file_path" in first):
            return False

        return True

    # --------- "Turn this config into domain objects" ---------

    def to_config_and_datasets(
        self, raw: Mapping[str, Any]
    ) -> Tuple[GlobalConfig, List[Dataset]]:
        """
        Map the JSON structure into:
        - GlobalConfig     (for UI / metadata)
        - list[Dataset]    (for the app to work with)
        """
        global_part = raw.get("global", {}) or {}
        dataset_entries = raw.get("datasets", []) or []

        # --- GlobalConfig (UI metadata only for now) ---
        global_config = GlobalConfig(
            ui_title=global_part.get("ui_title", "Single-cell Browser"),
            default_group=global_part.get("default_group", "Default"),
        )

        # --- Datasets (domain layer) ---
        datasets: List[Dataset] = []
        for i, entry in enumerate(dataset_entries):
            # Dataset.from_config_entry already knows how to handle
            # "file" vs "file_path" etc. – we just pass the dict through.
            ds = Dataset.from_config_entry(entry)
            datasets.append(ds)

        return global_config, datasets