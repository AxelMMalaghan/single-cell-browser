from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Dict, Type, Any

from sc_browser.core.dataset import Dataset

class BaseConfigAdapter(ABC):
    """
    Takes a raw config entry (from JSON) and returns a dict compatible with
    {@link Dataset.from_config_entry()}
    """
    id: str

    @abstractmethod
    def can_handle(self):
        """
        Returns true if this adapter can handle this config entry
        """
        raise NotImplementedError


    @abstractmethod
    def build_dataset(self, entry: Dict[str, Any]) -> Dataset:
        """
        Builds the dataset from the config entry
        """
        raise NotImplementedError
