from __future__ import annotations

import json
from pathlib import Path
from typing import List, Tuple, Set

from sc_browser.config.model import DatasetConfig, GlobalConfig
from sc_browser.core.dataset_loader import from_config
from sc_browser.core.dataset import Dataset


def load_global_config(path: Path) -> GlobalConfig:
    """
    Loads the top-level application configuration from a JSON file or directory

    This supports two layout styles:

    1 - Directory
    2 - Legacy (single file) - Only here for demo datasets TODO: migrate away after demo on Friday
    :param path: Points to a directory or single file containing global.json AND datasets.json
    :return: Parsed 'GlobalConfig' instance describing app setup
    :raises FileNotFoundError: if path or file does not exist
    """

    if path.is_dir():
        return _load_global_from_dir(path)
    if path.is_file():
        return _load_global_from_file(path)
    raise FileNotFoundError(f"File not found at {path}")


def load_datasets(path: Path) -> Tuple[GlobalConfig, List[Dataset]]:
    """
    Loads the global configuration and instantiates all Dataset objects

    Main entrypoint used by UI/services

    1. Loads the top-level 'GlobalConfig' from 'path' supporting both directory and file.
    2. Merges:
        - All 'DatasetConfig' instances defined in the config and
        - Any additional datasets discovered under 'GlobalConfig.data_root'
    3. Materialises each 'DatasetConfig' into a 'Dataset' by calling 'Dataset.from_config'
    :param path: Path to config directory
    :return: A tuple of (GlobalConfig, Datasets)
    """
    global_config = load_global_config(path)

    # Merge configured + discovered
    all_cfgs: List[DatasetConfig] = [
        *global_config.datasets,
        *_discover_dataset_configs(global_config),
    ]

    datasets: List[Dataset] = []
    for ds_cfg in all_cfgs:
        datasets.append(from_config(ds_cfg))

    return global_config, datasets



def _load_global_from_dir(root: Path) -> GlobalConfig:
    """
    Loads configuration from a directory using the new multi-file layout

    Expected structure:

        root/
            global.json
            datasets/
                dataset_1.json
                dataset_2.json

    Each file in 'datasets/' is parsed into a 'DatasetConfig'. The resulting 'Global Config' includes:

    - ui_title: title for UI, defaults to 'Single-Cell Browser'
    - default_group: group for UI, defaults to 'Default'
    - datasets: list of DatasetConfigs
    - data_root: root directory for datasets

    :param root: Directory containing 'global.json' and datasets/
    :return: A GlobalConfig instance
    :raises FileNotFoundError: global.json does not exist
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
            data_root=Path(raw_global.get("data_root", "data")).resolve()
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

def _discover_dataset_configs(global_config: GlobalConfig) -> List[DatasetConfig]:
    """
    Find .h5ad files under data_root that are not already configured
    and create minimal DatasetConfig entries for them.
    """
    if global_config.data_root is None:
        return []

    data_root = global_config.data_root
    if not data_root.is_dir():
        return []

    # Paths already covered by explicit configs
    configured_paths: Set[Path] = {
        cfg.path.resolve() for cfg in global_config.datasets
    }

    discovered: List[DatasetConfig] = []
    idx_offset = len(global_config.datasets)

    for idx, path in enumerate(sorted(data_root.rglob("*.h5ad"))):
        if path.resolve() in configured_paths:
            continue

        raw = {
            "name": path.stem,
            "group": "Discovered",
            "path": str(path),
            # minimal config – no obs_columns yet
            # views will still work for things that don't need cluster/condition
        }
        discovered.append(
            DatasetConfig.from_raw(raw, source_path=path, index=idx_offset + idx)
        )

    return discovered