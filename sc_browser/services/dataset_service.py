from __future__ import annotations

import logging
import os
from collections import OrderedDict
from typing import Dict, Iterator, Mapping

from sc_browser.core.configs import DatasetConfig
from sc_browser.core.dataset import Dataset
from sc_browser.config.dataset_loader import from_config, DatasetConfigError

logger = logging.getLogger(__name__)


class DatasetManager(Mapping[str, Dataset]):
    """
    Central service for managing datasets.
    Implements LRU caching to prevent unbounded memory usage.
    """

    # Default to keeping 5 large datasets in memory.
    # Can be tuned via env var SC_BROWSER_MAX_LOADED_DATASETS
    DEFAULT_MAX_LOADED = 7

    def __init__(self, cfg_by_name: Dict[str, DatasetConfig]):
        self._cfg_by_name = cfg_by_name
        # Use OrderedDict for LRU behavior: oldest items are at the beginning
        self._loaded: OrderedDict[str, Dataset] = OrderedDict()

        try:
            self._max_loaded = int(
                os.getenv(
                    "SC_BROWSER_MAX_LOADED_DATASETS", str(self.DEFAULT_MAX_LOADED)
                )
            )
        except ValueError:
            self._max_loaded = self.DEFAULT_MAX_LOADED

    def __getitem__(self, name: str) -> Dataset:
        # 1. Fast path: already loaded (LRU Cache Hit)
        if name in self._loaded:
            # Move to end to mark as recently used
            self._loaded.move_to_end(name)
            return self._loaded[name]

        # 2. Check config existence
        cfg = self._cfg_by_name.get(name)
        if cfg is None:
            raise KeyError(f"Unknown dataset '{name}'")

        # 3. Lazy load
        try:
            logger.info(
                "Lazy-loading dataset",
                extra={"dataset": cfg.name, "path": str(getattr(cfg, "path", ""))},
            )
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

        # 4. Cache Management (LRU Eviction)
        self._loaded[name] = ds

        if len(self._loaded) > self._max_loaded:
            # Pop the first item (FIFO/Oldest)
            evicted_name, _ = self._loaded.popitem(last=False)
            logger.info(
                "Evicted dataset from memory cache to free space",
                extra={"dataset": evicted_name},
            )

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
        Clears the loaded cache to ensure consistency.
        """
        self._cfg_by_name = new_cfg_by_name
        logger.info("Dataset config refreshed. Clearing memory cache.")
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

    def get(self, key: str, default=None) -> Dataset | None:
        try:
            return self[key]
        except KeyError:
            return default

    def refresh(self, name_by_key: Dict[str, str]) -> None:
        """
        Update the key-to-name mapping (e.g. after configuration reload).
        """
        self._name_by_key = name_by_key
