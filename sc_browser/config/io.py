from __future__ import annotations
from typing import List
from pathlib import Path

import json

from sc_browser.config.model import GlobalConfig, DatasetConfig

def load_global_config(path: Path) -> GlobalConfig:
    """
    Load and parse the top-level config file.

    Purpose:
    - Reads a JSON config file
    - Normalises each dataset entry into a {@link DatasetConfig}
    - Returns a {@link GlobalConfig} object that the rest of the app can work with

    Expected schema (high-level):
        {
          "ui_title": "Single-cell browser",
          "default_group": "some-group",
          "datasets": [
            {
              "name": "...",
              "group": "...",
              "path": "...",
              "obs": {
                "cluster": "...",
                "condition": "...",
                "sample": "..."
              }
            },
            ...
          ]
        }

    :param path: filesystem path to the JSON config file
    :return: a fully constructed {@link GlobalConfig} instance
    :raises json.JSONDecodeError: if the file is not valid JSON
    :raises OSError: if the file cannot be opened
    """
    with path.open() as f:
        raw_config = json.load(f)

    datasets: List[DatasetConfig] = []
    for idx, entry in enumerate(raw_config.get("datasets", [])):
        datasets.append(DatasetConfig.from_raw(entry, source_path=path, index=idx))

    return GlobalConfig(
        ui_title=raw_config.get("ui_title", "Single-cell browser"),
        default_group=raw_config.get("default_group", "group"),
        datasets=datasets
    )




def load_datasets(path: Path) -> List[DatasetConfig]:
    """
    Helper to load only the dataset config

    Purpose:
    - For callers that only care about the list of datasets and do not need the full {@link GlobalConfig}
    - Keeps call sites simple: they just pass a config path and get back a list of {@link DatasetConfig} instances

    Behaviour:
    - Delegates to {@link load_global_config(path)}
    - Returns the 'datasets' list from the resulting {@link GlobalConfig}

    :param path: filesystem path to the JSON config file
    :return: a list of {@link DatasetConfig} entries parsed from the config
    """
    return load_global_config(path).datasets