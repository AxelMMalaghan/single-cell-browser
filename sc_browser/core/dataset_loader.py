from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import List

from sc_browser.config.model import DatasetConfig, GlobalConfig

logger = logging.getLogger(__name__)


def load_global_config(config_root: Path) -> GlobalConfig:
    """
    Load the global UI and dataset configuration.
    """
    path = config_root / "global.json"
    logger.info("Loading global config", extra={"config_root": str(config_root)})

    if not path.is_file():
        # Minimal default if missing
        return GlobalConfig(
            ui_title="Single Cell Browser",
            default_group="Default",
            datasets=scan_datasets(config_root),
        )

    with open(path, "r") as f:
        raw = json.load(f)

    # If the global config doesn't list datasets, scan the directory
    datasets = raw.get("datasets")
    if datasets is None:
        datasets = scan_datasets(config_root)
    else:
        datasets = [
            DatasetConfig.from_raw(d, path, i) for i, d in enumerate(datasets)
        ]

    return GlobalConfig(
        ui_title=raw.get("ui_title", "Single Cell Browser"),
        default_group=raw.get("default_group", "Default"),
        datasets=datasets,
    )


def scan_datasets(config_root: Path) -> List[DatasetConfig]:
    """
    Scan the config/datasets directory for .json files and return a list of DatasetConfigs.

    Includes a fix to ignore macOS hidden metadata files (._*).
    """
    dataset_dir = config_root / "datasets"
    logger.info("Scanning for dataset configurations in: %s", dataset_dir)

    configs = []
    if not dataset_dir.is_dir():
        return configs

    for i, p in enumerate(sorted(dataset_dir.glob("*.json"))):
        # FIX: Ignore macOS hidden 'Apple Double' files which cause decoding errors
        if p.name.startswith("._"):
            continue

        logger.info("Loading dataset config: %s", p.name)
        try:
            with open(p, "r") as f:
                raw = json.load(f)
            configs.append(DatasetConfig.from_raw(raw, p, i))
        except Exception as e:
            logger.error("Failed to load %s: %s", p.name, e)

    return configs