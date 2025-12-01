from __future__ import annotations
from typing import List
from pathlib import Path

import json

from .config_model import GlobalConfig, DatasetConfig
from .core.dataset import Dataset
from .importing.adapter_registry import create_single_cell_registry


def load_datasets(config_path: str | Path) -> tuple[GlobalConfig, List[Dataset]]:
    """
    Load dataset config and instantiate Dataset objects.
    """
    config_path = Path(config_path)
    raw = json.loads(config_path.read_text())

    global_part = raw.get("global", {})
    dataset_entries = raw.get("datasets", [])

    global_config = GlobalConfig(
        ui_title=global_part.get("ui_title", "Single-Cell Browser"),
        default_group=global_part.get("default_group", "Default"),
        datasets=[
            DatasetConfig(raw=entry, source_path=config_path, index=i)
            for i, entry in enumerate(dataset_entries)
        ],
    )

    registry = create_single_cell_registry()

    datasets: List[Dataset] = [
        registry.create_dataset(ds_cfg.raw)
        for ds_cfg in global_config.datasets
    ]

    return global_config, datasets