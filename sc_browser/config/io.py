from __future__ import annotations

import json
from pathlib import Path
from typing import List

from sc_browser.config import GlobalConfig
from .model import GlobalConfig, DatasetConfig
from ..core import Dataset


def load_global_config(path: Path) -> GlobalConfig:
    """
    Loads the top level application configuration from a JSON file
    :param path:
    :return:
    """

    if path.is_dir():
        return _load_global_from_dir(path)
    if path.is_file():
        return _load_global_from_file(path)
    raise FileNotFoundError(f"File not found at {path}")

def load_datasets(path: Path) -> List[DatasetConfig]:
    """
    Load dataset configs from the given config path
    :param path: can either be a legacy JSON file or a config root directory
    :return: loaded config datasets
    """
    global_config = load_global_config(path)

    datasets: List[Dataset] = []
    for ds_cfg in global_config.datasets:
        datasets.append(
            Dataset.from_config(ds_cfg)
        )

    return global_config, datasets





def _load_global_from_dir(root: Path) -> GlobalConfig:
    """
    New loader

    :param root: is a directory containing :
                    root/global.json
                    root/dataset/*.json
    :return: the Global config and a list of Datasets
    """

    global_path = root / "global.json"
    if not global_path.is_file():
        raise FileNotFoundError(f"File not found at {global_path}")

    with global_path.open() as f:
        raw_global = json.load(f)

    datasets_dir = root / "datasets"
    datasets: List[DatasetConfig] = []

    if datasets_dir.is_dir():
        for idx, config_file in enumerate(sorted(datasets_dir.glob("*.json"))):
            with config_file.open() as f:
                raw = json.load(f)
            datasets.append(
                DatasetConfig.from_raw(raw, source_path=config_file, index=idx)
            )
        return GlobalConfig(
            ui_title=raw_global.get("ui_title", "Single-Cell Browser"),
            default_group=raw_global.get("default_group", "Default"),
            datasets=datasets,
        )


def _load_global_from_file(path: Path) -> GlobalConfig:
    """
    Legacy loader: single JSON file containing both global and datasets.

    {
      "ui_title": "...",
      "default_group": "...",
      "datasets": [ {...}, {...} ]
    }
    """
    with path.open() as f:
        raw_config = json.load(f)

    datasets: List[DatasetConfig] = []
    for idx, entry in enumerate(raw_config.get("datasets", [])):
        datasets.append(DatasetConfig.from_raw(entry, source_path=path, index=idx))

    return GlobalConfig(
        ui_title=raw_config.get("ui_title", "Single-cell browser"),
        default_group=raw_config.get("default_group", "group"),
        datasets=datasets,
    )
