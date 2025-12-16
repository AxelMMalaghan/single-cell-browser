from __future__ import annotations

import logging
from typing import Dict, Iterator, Mapping

from sc_browser.config.model import DatasetConfig
from sc_browser.core.dataset import Dataset
from sc_browser.core.dataset_loader import from_config, DatasetConfigError

logger = logging.getLogger(__name__)


class DatasetManager(Mapping[str, Dataset]):
    """
    Central service for managing datasets.
    Implements the Mapping interface (dict-like) to support lazy loading
    transparency for the UI layer.
    """

    def __init__(self, cfg_by_name: Dict[str, DatasetConfig]):
        self._cfg_by_name = cfg_by_name
        self._loaded: Dict[str, Dataset] = {}

    def __getitem__(self, name: str) -> Dataset:
        # 1. Fast path: already materialised
        if name in self._loaded:
            return self._loaded[name]

        # 2. Check config existence
        cfg = self._cfg_by_name.get(name)
        if cfg is None:
            raise KeyError(f"Unknown dataset '{name}'")

        # 3. Lazy load
        try:
            logger.info("Lazy-loading dataset", extra={"dataset": cfg.name, "path": str(getattr(cfg, "path", ""))})
            ds = from_config(cfg)
        except DatasetConfigError as e:
            logger.error(
                "Dataset config error on load",
                extra={"dataset": cfg.name, "error": str(e)},
            )
            raise
        except Exception:
            logger.exception(
                "Unexpected error while loading dataset",
                extra={"dataset": cfg.name},
            )
            raise

        self._loaded[name] = ds
        return ds

    def __iter__(self) -> Iterator[str]:
        return iter(self._cfg_by_name)

    def __len__(self) -> int:
        return len(self._cfg_by_name)

    def get(self, name: str, default=None) -> Dataset | None:
        try:
            return self[name]
        except KeyError:
            return default

    def is_loaded(self, name: str) -> bool:
        return name in self._loaded

    def refresh_config(self, new_cfg_by_name: Dict[str, DatasetConfig]) -> None:
        """
        Update the configuration map (e.g. after import).
        Clears the loaded cache for any datasets that have changed or been removed.
        """
        self._cfg_by_name = new_cfg_by_name
        # Conservative: clear everything to ensure consistency
        # Optimization: could diff keys, but safely reloading is better
        self._loaded.clear()


class DatasetKeyManager(Mapping[str, Dataset]):
    """
    Lazy view of datasets by their stable 'key'.
    Delegates actual loading to the parent DatasetManager.
    """

    def __init__(self, manager: DatasetManager, name_by_key: Dict[str, str]):
        self._manager = manager
        self._name_by_key = name_by_key

    def __getitem__(self, key: str) -> Dataset:
        if key not in self._name_by_key:
            raise KeyError(f"Unknown dataset key '{key}'")
        name = self._name_by_key[key]
        return self._manager[name]

    def __iter__(self) -> Iterator[str]:
        return iter(self._name_by_key)

    def __len__(self) -> int:
        return len(self._name_by_key)