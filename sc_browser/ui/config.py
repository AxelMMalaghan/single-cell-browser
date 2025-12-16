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

    # Lazy registry / config listing (no materialisation)
    dataset_names: List[str] = field(default_factory=list)

    # Materialised datasets ONLY
    datasets: List[Dataset] = field(default_factory=list)
    dataset_by_name: Dict[str, Dataset] = field(default_factory=dict)
    dataset_by_key: Dict[str, Dataset] = field(default_factory=dict)
    default_dataset: Optional[Dataset] = None

    registry: ViewRegistry = None  # or keep your existing init pattern
    export_service: ExportService = None
    enable_dataset_management: bool = False
