from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import List, Tuple

from sc_browser.config.model import DatasetConfig, GlobalConfig
from sc_browser.core.dataset_loader import from_config, DatasetConfigError
from sc_browser.core.dataset import Dataset

logger = logging.getLogger(__name__)


def load_global_config(root: Path) -> GlobalConfig:
    """
    Load configuration from a directory using the multi-file layout.

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
                 (Auto-discovery may be handled elsewhere; this function just
                  parses what is explicitly configured here.)

    :param root: Directory containing 'global.json' and optionally 'datasets/'.
    :return: A GlobalConfig instance.
    :raises FileNotFoundError: if global.json does not exist.
    """
    logger.info(
        "Loading global config",
        extra={"config_root": str(root)},
    )

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


def load_datasets(path: Path) -> Tuple[GlobalConfig, List[Dataset]]:
    """
    Load the global configuration and instantiate all Dataset objects.

    Main entrypoint used by UI/services.

    1. Loads the top-level GlobalConfig from 'path'.
    2. Iterates over GlobalConfig.datasets and calls `from_config` for each.
    3. Skips any dataset whose config is invalid, logging the error.
    4. Returns only successfully materialised Dataset objects.

    :param path: Path to config directory.
    :return: A tuple of (GlobalConfig, List[Dataset]).
    :raises RuntimeError: if no valid datasets could be loaded.
    """
    global_config = load_global_config(path)

    all_cfgs: List[DatasetConfig] = global_config.datasets

    datasets: List[Dataset] = []
    failed = 0

    for ds_cfg in all_cfgs:
        try:
            ds = from_config(ds_cfg)
        except DatasetConfigError as e:
            failed += 1
            logger.error(
                "Skipping dataset due to config error",
                extra={
                    "dataset": ds_cfg.name,
                    "path": str(getattr(ds_cfg, "path", "")),
                    "error": str(e),
                },
            )
            continue
        except Exception as e:  # hard guard against unexpected issues
            failed += 1
            logger.exception(
                "Unexpected error while loading dataset; skipping",
                extra={
                    "dataset": ds_cfg.name,
                    "path": str(getattr(ds_cfg, "path", "")),
                },
            )
            continue

        datasets.append(ds)

    logger.info(
        "Datasets loaded from config root",
        extra={
            "config_root": str(path),
            "n_datasets": len(datasets),
            "n_failed": failed,
            "dataset_names": [ds.name for ds in datasets],
        },
    )

    if not datasets:
        # At this point, the app would be useless anyway, so fail fast with a clear error.
        raise RuntimeError(
            f"No valid datasets could be loaded from config root: {path}"
        )

    return global_config, datasets