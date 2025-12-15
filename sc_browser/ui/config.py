from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

from sc_browser.config.model import GlobalConfig
from sc_browser.core.dataset import Dataset
from sc_browser.core.view_registry import ViewRegistry
from sc_browser.services.metadata_export_service import ExportService


@dataclass
class AppConfig:
    """
    Holds shared state for the Dash app: config root, loaded datasets and
    the view registry. This is passed into layout + callback registration
    functions instead of using module-level globals.
    """
    config_root: Path
    global_config: GlobalConfig
    datasets: List[Dataset]
    dataset_by_name: Dict[str, Dataset]
    dataset_by_key: Dict[str, Dataset]
    default_dataset: Optional[Dataset]
    registry: ViewRegistry
    export_service: ExportService
    enable_dataset_management: bool = False  # default safe for prod