from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import List, Dict

from sc_browser.config.model import DatasetConfig, GlobalConfig

logger = logging.getLogger(__name__)


def load_global_config(root: Path) -> GlobalConfig:
    """
    Load configuration from a directory using the multi-file layout.
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
        # --- ADDED LOGGING HERE ---
        logger.info(f"Scanning for dataset configurations in: {datasets_dir}")

        files = sorted(list(datasets_dir.glob("*.json")))
        if not files:
            logger.warning(f"No .json files found in {datasets_dir}")

        for idx, config_file in enumerate(files):
            logger.info(f"Loading dataset config: {config_file.name}")
            try:
                with config_file.open() as f:
                    raw = json.load(f)
                datasets.append(
                    DatasetConfig.from_raw(raw, source_path=config_file, index=idx)
                )
            except Exception as e:
                logger.error(f"Failed to load {config_file.name}: {e}")
        # --------------------------
    else:
        logger.warning(f"Datasets directory not found at: {datasets_dir}")

    # Resolve data_root properly:
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


def load_dataset_registry(path: Path) -> tuple[GlobalConfig, Dict[str, DatasetConfig]]:
    """
    Load global config + dataset config objects only (NO AnnData loading).
    Returns mapping of dataset name -> DatasetConfig.
    """
    global_config = load_global_config(path)
    all_cfgs: List[DatasetConfig] = global_config.datasets

    cfg_by_name: Dict[str, DatasetConfig] = {}
    duplicates: List[str] = []

    for ds_cfg in all_cfgs:
        if ds_cfg.name in cfg_by_name:
            duplicates.append(ds_cfg.name)
            continue
        cfg_by_name[ds_cfg.name] = ds_cfg

    if duplicates:
        raise RuntimeError(f"Duplicate dataset names in config: {sorted(set(duplicates))}")

    if not cfg_by_name:
        # Changed from RuntimeError to warning so app can start empty if needed
        logger.warning(f"No datasets configured under: {path}")

    logger.info(
        "Dataset registry loaded (lazy mode; datasets not materialised)",
        extra={
            "config_root": str(path),
            "n_dataset_configs": len(cfg_by_name),
            "dataset_names": sorted(cfg_by_name.keys()),
        },
    )

    return global_config,
