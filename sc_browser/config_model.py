from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List


@dataclass
class DatasetConfig:
    """
    C
    """
    raw: Dict[str, Any]
    source_path: Path
    index: int

    @property
    def name(self) -> str:
        return self.raw.get("name", f"Dataset {self.index}")


@dataclass
class GlobalConfig:
    ui_title: str
    default_group: str
    datasets: List[DatasetConfig]