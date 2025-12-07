from __future__ import annotations

from typing import List, Any, Dict

from sc_browser.core.dataset import Dataset
from sc_browser.core.base_adapter import BaseConfigAdapter
from sc_browser.import_adapters.mapped_anndata_adapter import MappedAnnDataAdapter


class AdapterRegistry:
    """
    Registry of adapters for turning dataset config entries into Dataset objects.

    Purpose:
    - Acts as the source of truth for which adapters are registered
    - Decouples io/config parsing from adapter implementations by exposing {@link create_dataset_(entry) as the construction API}
    - Enables support for multiple config schemas or data sources without changing call sites

    Design Notes:
    - Stores instances of {@link BaseConfigAdapter}
    - Enforces variants:
        * only {@link BaseConfigAdapter} instances are registered
        * each adapter's 'id' is unique
    - Resolution is first-hit: the first adapters where {@link BaseConfigAdapter} can handle the entry
    - Centralising creation logic makes it easier to extend later without touching rest of codebase
    """

    def __init__(self) -> None:
        self._adapters: List[BaseConfigAdapter] = []

    @property
    def adapters(self) -> List[BaseConfigAdapter]:
        """
        :return: a list of all adapters
        """
        return self._adapters

    def register(self, adapter: BaseConfigAdapter) -> None:
        """
        Register a new adapter instance.
        :raises ValueError: if an adapter with the same id already exists.
        """
        for existing in self._adapters:
            if existing.id == adapter.id:
                raise ValueError(f"The adapter '{existing.id}' already exists")

        self._adapters.append(adapter)

    def create_dataset(self, entry: Dict[str, Any]) -> Dataset:
        """
        Finds the adapter that can handle the given entry
        and uses it to build a {@link Dataset} object

        :param entry: the incoming dataset config dict
        :return: a {@link Dataset} object
        :raises ValueError: if no adapter can handle the entry
        """
        for adapter in self._adapters:
            if adapter.can_handle(entry):
                return adapter.build_dataset(entry)

        name = entry.get("name", "<unnamed dataset>")
        schema = entry.get("schema", "<no schema specified>")
        raise ValueError(
            f"No config adapter could handle dataset entry "
            f"name='{name}', schema='{schema}'. "
            f"Registered adapters: {[a.id for a in self._adapters]}"
        )

    def get_adapters(self) -> List[BaseConfigAdapter]:
        """
        :return: a list of all adapters (shallow copy)
        """
        return list(self._adapters)


def create_single_cell_registry() -> AdapterRegistry:
    """
    Builds a registry with all known single-cell adapters.
    For now, just the simple entry-based adapter.
    """
    registry = AdapterRegistry()
    registry.register(MappedAnnDataAdapter())
    return registry

