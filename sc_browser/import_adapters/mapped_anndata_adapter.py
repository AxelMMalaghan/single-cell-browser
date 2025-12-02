from pathlib import Path
from typing import Any, Dict

from sc_browser.core.dataset import Dataset
from sc_browser.core.base_adapter import BaseConfigAdapter
from sc_browser.config.model import DatasetConfig
from sc_browser.core.adapters.anndata_adapter import AnnDataAdapter

class MappedAnnDataAdapter(BaseConfigAdapter):

    id = "anndata_mapped"

    def can_handle(self, entry: Dict[str, Any]) -> bool:
        return entry.get("schema") == self.id and "file" in entry and "obs_columns" in entry

    def build_dataset(self, entry: Dict[str, Any]) -> Dataset:
        config = DatasetConfig(
            raw=entry,
            source_path=Path(entry.get("file", ".")),
            index=0  # runtime value, can be updated if used in batch configs
        )
        adapter = AnnDataAdapter(config)
        return Dataset(
            name=config.name,
            adapter=adapter,
            config=config
        )
