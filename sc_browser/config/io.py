from __future__ import annotations

import json
from pathlib import Path
from typing import List, Tuple

from .model import GlobalConfig, DatasetConfig
from sc_browser.core.dataset import Dataset


def load_global_config(path: str | Path) -> GlobalConfig:
    path = Path(path)
    raw = json.loads(path.read_text())

    global_part = raw.get("global", {}) or {}
    dataset_entries = raw.get("datasets", []) or []

    datasets_cfg: List[DatasetConfig] = [
        DatasetConfig(raw=entry, source_path=path, index=i)
        for i, entry in enumerate(dataset_entries)
    ]

    return GlobalConfig(
        ui_title=global_part.get("ui_title", "Single-Cell Browser"),
        default_group=global_part.get("default_group", "Default"),
        datasets=datasets_cfg,
    )


def load_datasets(path: str | Path) -> Tuple[GlobalConfig, List[Dataset]]:
    global_config = load_global_config(path)

    datasets: List[Dataset] = [
        Dataset.from_config_entry(ds_cfg.raw)
        for ds_cfg in global_config.datasets
    ]

    return global_config, datasets