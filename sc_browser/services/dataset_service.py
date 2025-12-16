from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict

from sc_browser.config.model import DatasetConfig
from sc_browser.core.dataset import Dataset
from sc_browser.core.dataset_loader import from_config, DatasetConfigError

logger = logging.getLogger(__name__)


@dataclass
class DatasetManager:
    cfg_by_name: Dict[str, DatasetConfig]
    _loaded: Dict[str, Dataset] = field(default_factory=dict)

    def get(self, name: str) -> Dataset:
        # Fast path
        ds = self._loaded.get(name)
        if ds is not None:
            return ds

        cfg = self.cfg_by_name.get(name)
        if cfg is None:
            raise KeyError(f"Unknown dataset '{name}'")

        try:
            logger.info("Lazy-loading dataset", extra={"dataset": cfg.name, "path": str(getattr(cfg, "path", ""))})
            ds = from_config(cfg)
        except DatasetConfigError as e:
            # In lazy mode, we *shouldn't* skip silently — user explicitly asked for this dataset.
            logger.error(
                "Dataset config error on load",
                extra={"dataset": cfg.name, "path": str(getattr(cfg, "path", "")), "error": str(e)},
            )
            raise
        except Exception:
            logger.exception(
                "Unexpected error while loading dataset",
                extra={"dataset": cfg.name, "path": str(getattr(cfg, "path", ""))},
            )
            raise

        self._loaded[name] = ds
        return ds

    def is_loaded(self, name: str) -> bool:
        return name in self._loaded

    def loaded_names(self) -> list[str]:
        return sorted(self._loaded.keys())
