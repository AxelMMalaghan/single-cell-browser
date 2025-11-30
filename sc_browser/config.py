from __future__ import annotations
import json
from pathlib import Path
from typing import List, Tuple

from .config_model import DatasetConfig, GlobalConfig
from .core.dataset import Dataset
from .importing.adapter_registry import AdapterRegistry


def load_datasets(config_path: str | Path):
    """
    Load dataset config and instantiate Dataset objects
    """

    raw = json.loads(Path(config_path).read_text())

    global_part = raw.get("global", {})
    dataset_entries = raw.get("datasets", [])

    global_config = GlobalConfig(
        ui_title=global_part.get("ui_title", "Single-Cell Browser"),
        default_group=global_part.get("default_group", "Default"),
        datasets = [
            DatasetConfig(raw=entry, source_path=Path(config_path), index=i)
            for i, entry in enumerate(dataset_entries)
        ],
    )

    registry = AdapterRegistry.default()

    datasets = []
    for ds_config in global_config.datasets:
        adapter = registry.get_adapter(ds_config.raw)
        dataset = adapter.choose_adapter(ds_config.raw)
        datasets.append(dataset)

    return global_config, datasets



