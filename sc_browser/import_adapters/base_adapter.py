from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Dict, Any

from sc_browser.core.dataset import Dataset

class BaseConfigAdapter(ABC):
    """
    Abstract Base Class for all config adapters

    Defines the contract that every adapter in the app must follow
    - expose an 'id' used internally to identify the adapter
    - implement 'can_handle' method to check if this adapter is able to standardise this config
    - implement 'build_dataset' method to build the dataset from the config entry
    """
    id: str

    @abstractmethod
    def can_handle(self, entry: Dict[str, Any]) -> bool:
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



    # second iteration, converts .h5ad files


    # third iteration - legacy data formats (seurat) converts, etc...

