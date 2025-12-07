from __future__ import annotations

import json
from pathlib import Path
from typing import List, Tuple, Set

from sc_browser.config.model import DatasetConfig, GlobalConfig
from sc_browser.core.dataset_loader import from_config
from sc_browser.core.dataset import Dataset


def load_global_config(path: Path) -> GlobalConfig:
    """
    Load the top-level application configuration from a JSON file or directory.

    Supports two layouts:

    1. Directory layout (preferred):
        root/
            global.json
            datasets/
                dataset_1.json
                dataset_2.json
                ...

    2. Legacy single-file layout (for demo / backwards compatibility):
        config.json
        {
          "ui_title": "...",
          "default_group": "...",
          "datasets": [ {...}, {...} ]
        }

    :param path: Directory containing global.json + datasets/ OR a single legacy JSON file.
    :return: Parsed GlobalConfig instance describing app setup.
    :raises FileNotFoundError: if path or file does not exist.
    """
    if path.is_dir():
        return _load_global_from_dir(path)
    if path.is_file():
        return _load_global_from_file(path)
    raise FileNotFoundError(f"File or directory not found at {path}")


def load_datasets(path: Path) -> Tuple[GlobalConfig, List[Dataset]]:
    """
    Load the global configuration and instantiate all Dataset objects.

    Main entrypoint used by UI/services.

    1. Loads the top-level GlobalConfig from 'path' (directory or file).
    2. Merges:
        - All DatasetConfig instances defined explicitly in the config, and
        - Any additional datasets discovered under GlobalConfig.data_root.
    3. Materialises each DatasetConfig into a Dataset by calling 'dataset_loader.from_config'.

    :param path: Path to config directory or legacy config file.
    :return: A tuple of (GlobalConfig, List[Dataset]).
    """
    global_config = load_global_config(path)

    # Merge configured + discovered configs (discovery is a no-op if data_root is None)
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
    Load configuration from a directory using the new multi-file layout.

    Expected structure:

        root/
            global.json
            datasets/
                dataset_1.json
                dataset_2.json
                ...

    Each file in 'datasets/' is parsed into a DatasetConfig. The resulting GlobalConfig includes:

    - ui_title: title for UI, defaults to 'Single-Cell Browser'
    - default_group: group for UI, defaults to 'Default'
    - datasets: list of DatasetConfigs
    - data_root: root directory for datasets, used for auto-discovery of .h5ad files.
                 If 'data_root' is relative in global.json, it is resolved relative to 'root'.

    :param root: Directory containing 'global.json' and optionally 'datasets/'.
    :return: A GlobalConfig instance.
    :raises FileNotFoundError: if global.json does not exist.
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

    # Resolve data_root properly:
    # - Absolute paths are used as-is.
    # - Relative paths are resolved relative to the config root directory.
    data_root_raw = raw_global.get("data_root")
    if data_root_raw is None:
        data_root = None
    else:
        data_root_path = Path(data_root_raw)
        if data_root_path.is_absolute():
            data_root = data_root_path
        else:
            data_root = (root / data_root_path).resolve()

    return GlobalConfig(
        ui_title=raw_global.get("ui_title", "Single-Cell Browser"),
        default_group=raw_global.get("default_group", "Default"),
        datasets=datasets,
        data_root=data_root,
    )


def _load_global_from_file(path: Path) -> GlobalConfig:
    """
    Legacy loader: single JSON file containing both global and dataset configuration.

    Example:

    {
      "ui_title": "...",
      "default_group": "...",
      "datasets": [ {...}, {...} ]
    }

    NOTE:
    - Legacy layout does not set 'data_root' explicitly; GlobalConfig should default it to None.
    - Dataset paths in this mode are interpreted by DatasetConfig.from_raw.
    """
    with path.open() as f:
        raw_config = json.load(f)

    datasets: List[DatasetConfig] = []

    for idx, entry in enumerate(raw_config.get("datasets", [])):
        # Minimal sanity check; tighten if you want stricter validation.
        if "name" not in entry:
            raise ValueError(f"Dataset config at index {idx} is missing required field 'name'")
        datasets.append(DatasetConfig.from_raw(entry, source_path=path, index=idx))

    return GlobalConfig(
        ui_title=raw_config.get("ui_title", "Single-Cell Browser"),
        default_group=raw_config.get("default_group", "group"),
        datasets=datasets,
        # data_root intentionally left to GlobalConfig default (likely None) in legacy mode
    )


def _discover_dataset_configs(global_config: GlobalConfig) -> List[DatasetConfig]:
    """
    Find .h5ad files under data_root that are not already configured,
    and create minimal DatasetConfig entries for them.

    Discovery is intentionally minimal:
    - name: derived from filename stem
    - group: 'Discovered'
    - path: full path to the file

    Views that require specific obs-columns (cluster/condition) may not work
    fully with these until they are explicitly configured.
    """
    if global_config.data_root is None:
        return []

    data_root = global_config.data_root
    if not data_root.is_dir():
        return []

    # Paths already covered by explicit configs
    configured_paths: Set[Path] = {
        cfg.path.resolve() for cfg in global_config.datasets
        if getattr(cfg, "path", None) is not None
    }

    discovered: List[DatasetConfig] = []
    idx_offset = len(global_config.datasets)

    for idx, path in enumerate(sorted(data_root.rglob("*.h5ad"))):
        resolved = path.resolve()
        if resolved in configured_paths:
            continue

        raw = {
            "name": path.stem,
            "group": "Discovered",
            "path": str(resolved),
            # Minimal config – no obs_columns yet.
            # Cluster/condition-based views must be defensive for these.
        }
        discovered.append(
            DatasetConfig.from_raw(raw, source_path=path, index=idx_offset + idx)
        )

    return discovered
