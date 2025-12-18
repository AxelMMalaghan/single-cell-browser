from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

from sc_browser.config.model import GlobalConfig
from sc_browser.core.dataset import Dataset
from sc_browser.core.view_registry import ViewRegistry
from sc_browser.services.metadata_export_service import ExportService


@dataclass
class AppConfig:
    config_root: Path
    global_config: GlobalConfig
    dataset_names: List[str] = field(default_factory=list)
    datasets: List[Dataset] = field(default_factory=list)
    dataset_by_name: Dict[str, Dataset] = field(default_factory=dict)
    dataset_by_key: Dict[str, Dataset] = field(default_factory=dict)
    default_dataset: Optional[Dataset] = None

    registry: Optional[ViewRegistry] = None
    export_service: Optional[ExportService] = None
    enable_dataset_management: bool = False

    def validate(self) -> None:
        """Ensure all required services are attached before the app starts."""
        if self.registry is None:
            raise RuntimeError("AppConfig.registry must be initialized.")
        if self.export_service is None:
            raise RuntimeError("AppConfig.export_service must be initialized.")