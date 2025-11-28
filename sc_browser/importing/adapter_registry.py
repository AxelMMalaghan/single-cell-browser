from __future__ import annotations

from typing import List, Any, Dict

from sc_browser.core import Dataset
from sc_browser.importing.adapters import BaseConfigAdapter


class AdapterRegistry:
    """
    Handles the registry of adapters.

    Given a config entry, it finds the adapter which returns True (can handle the schema)
    """

    def __init__(self) -> None:
        self.adapters: List[BaseConfigAdapter] = []


    def register(self, adapter: BaseConfigAdapter) -> None:
        """
        Register a new adapter instance
        :param adapter: the adapter to instantiate
        :raises ValueError: if an adapter with the same id already exists
        """
        for existing in self.adapters:
            if existing.id == adapter.id:
                raise ValueError(f"The adapter {existing} already exists")

        self.adapters.append(adapter)

    def adapters(self) -> List[BaseConfigAdapter]:
        """
        :return: a list of all adapters
        """
        return list(self.adapters)

    def create_dataset(self, entry: Dict[str, Any]) -> Dataset:
        """
        Finds the first adapter that can handle the given entry

        creates the dataset from an entry
        :param entry:
        :return:
        :raises ValueError: if not adapter can handle the entry
        """

        for adapter in self.adapters:
            if adapter.can_handle(entry):
                return adapter.build_dataset(entry)


        name=entry.get["name", "<unnamed dataset>"]
        schema=entry.get["shema", "<no schema specified>"]
        raise ValueError(
            f"No config adapter could handle dataset entry "
            f"name='{name}', schema='{schema}'. "
            f"Registered adapters: {[a.id for a in self._adapters]}"
        )


def create_single_cell_registry() -> AdapterRegistry:
    """
    Builds a registry with all known single-cell adapters
    """

    registry = AdapterRegistry()
    # registry.register() //TODO: add registries based on real data
    return registry