# sc_browser/config.py

from __future__ import annotations
import json
from pathlib import Path

from .config_model import GlobalConfig
from .core.dataset import Dataset
from .importing.adapter_registry import AdapterRegistry


def load_datasets(config_path: str | Path) -> tuple[GlobalConfig, list[Dataset]]:
    """
    Single entrypoint for the app:
    - Reads raw JSON
    - Uses the adapter registry
    - Returns a GlobalConfig + list[Dataset]
    """
    config_path = Path(config_path)
    raw = json.loads(config_path.read_text())

    registry = AdapterRegistry.default()
    adapter = registry.choose_adapter(raw)

    # Adapter is responsible for:
    #   - building GlobalConfig
    #   - building Dataset objects
    global_config, datasets = adapter.to_config_and_datasets(raw)

    return global_config, datasets