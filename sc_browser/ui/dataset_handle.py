from dataclasses import dataclass
from typing import Optional

from sc_browser.core.dataset import Dataset

@dataclass(frozen=True)
class DatasetHandle:
    name: str
    key: Optional[str] = None  # if you have stable keys
    _dataset: Optional[Dataset] = None

    def is_materialised(self) -> bool:
        return self._dataset is not None