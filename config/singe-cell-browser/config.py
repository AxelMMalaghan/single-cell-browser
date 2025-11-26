from __future__ import annotations
import json
from pathlib import Path
from typing import List, Tuple

from .core.dataset import Dataset

def load_datasets(config_path: str = "config/datasets.json") -> Tuple[dict, list[Dataset, ...]]:
    """
    Load dataset config and instantiate Dataset objects
    """

    raw = json.loads(Path(config_path).read_text())

    global_config = raw.get("config", {})
    dataset_entries = raw["datasets"]

    datasets: List[Dataset] = [
        Dataset.from_config_entry(entry) for entry in dataset_entries
    ]
    return global_config, datasets