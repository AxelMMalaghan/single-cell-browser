from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import List, Dict, Tuple

from sc_browser.core.configs import DatasetConfig, GlobalConfig

logger = logging.getLogger(__name__)


def load_global_config(root: Path) -> GlobalConfig:
    """
    Load configuration from a directory using the multi-file layout.
    """
    logger.info("Loading global config", extra={"config_root": str(root)})

    global_path = root / "global.json"
    if not global_path.is_file():
        # Fallback to defaults if global.json is missing
        raw_global = {}
    else:
        with global_path.open() as f:
            raw_global = json.load(f)

    datasets_dir = root / "datasets"
    datasets: List[DatasetConfig] = []

    if datasets_dir.is_dir():
        logger.info(f"Scanning for dataset configurations in: {datasets_dir}")
        # Sort files for deterministic loading
        files = sorted(list(datasets_dir.glob("*.json")))

        for idx, config_file in enumerate(files):
            # FIX: Ignore macOS hidden 'Apple Double' files (._*) to prevent decoding errors
            if config_file.name.startswith("._"):
                continue

            logger.info(f"Loading dataset config: {config_file.name}")
            try:
                with config_file.open() as f:
                    raw = json.load(f)
                datasets.append(
                    DatasetConfig.from_raw(raw, source_path=config_file, index=idx)
                )
            except Exception as e:
                logger.error(f"Failed to load {config_file.name}: {e}")
    else:
        logger.warning(f"Datasets directory not found at: {datasets_dir}")

    # Resolve data_root properly
    data_root_raw = raw_global.get("data_root")
    data_root = Path(data_root_raw) if data_root_raw else None
    if data_root and not data_root.is_absolute():
        data_root = (root / data_root).resolve()

    return GlobalConfig(
        ui_title=raw_global.get("ui_title", "Single-Cell Browser"),
        default_group=raw_global.get("default_group", "Default"),
        datasets=datasets,
        data_root=data_root,
    )


def load_dataset_registry(path: Path) -> Tuple[GlobalConfig, Dict[str, DatasetConfig]]:
    """
    Load global config + dataset config mapping.
    Expected by sc_browser/ui/dash_app.py
    """
    global_config = load_global_config(path)
    all_cfgs: List[DatasetConfig] = global_config.datasets

    cfg_by_name: Dict[str, DatasetConfig] = {}
    for ds_cfg in all_cfgs:
        if ds_cfg.name in cfg_by_name:
            logger.warning(f"Duplicate dataset name ignored: {ds_cfg.name}")
            continue
        cfg_by_name[ds_cfg.name] = ds_cfg

    # CRITICAL: Returns exactly two items for the unpack in dash_app.py
    return global_config, cfg_by_name