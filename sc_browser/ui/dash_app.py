from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Mapping, Iterator

import dash_bootstrap_components as dbc
from dash import Dash
import logging

from .config import AppConfig
from sc_browser.config.dataset_loader import load_dataset_registry
from sc_browser.core.view_registry import ViewRegistry
from sc_browser.core.dataset import Dataset
from sc_browser.services.metadata_export_service import ExportService
from sc_browser.ui.layout.build_layout import build_layout
from sc_browser.ui.callbacks.callbacks_explore import register_explore_callbacks
from sc_browser.ui.callbacks.callbacks_dataset_import import register_dataset_import_callbacks
from sc_browser.ui.callbacks.callbacks_reports import register_reports_callbacks
from sc_browser.ui.callbacks.callbacks_dataset_preview import register_dataset_preview_callbacks
from ..config.model import DatasetConfig
from ..core.dataset_loader import DatasetConfigError, from_config

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class _LazyDatasetsByName(Mapping[str, Dataset]):
    cfg_by_name: Dict[str, DatasetConfig]
    cache: Dict[str, Dataset]

    def __getitem__(self, name: str) -> Dataset:
        if name in self.cache:
            return self.cache[name]

        cfg = self.cfg_by_name[name]  # KeyError if unknown (good)
        try:
            logger.info("Lazy-loading dataset", extra={"dataset": cfg.name, "path": str(getattr(cfg, "path", ""))})
            ds = from_config(cfg)
        except DatasetConfigError:
            logger.exception("Dataset config error while lazy-loading", extra={"dataset": cfg.name})
            raise
        except Exception:
            logger.exception("Unexpected error while lazy-loading dataset", extra={"dataset": cfg.name})
            raise

        self.cache[name] = ds
        return ds

    def __iter__(self) -> Iterator[str]:
        return iter(self.cfg_by_name)

    def __len__(self) -> int:
        return len(self.cfg_by_name)


@dataclass(frozen=True)
class _LazyDatasetsByKey(Mapping[str, Dataset]):
    name_by_key: Dict[str, str]
    by_name: _LazyDatasetsByName

    def __getitem__(self, key: str) -> Dataset:
        name = self.name_by_key[key]
        return self.by_name[name]

    def __iter__(self) -> Iterator[str]:
        return iter(self.name_by_key)

    def __len__(self) -> int:
        return len(self.name_by_key)


def _build_view_registry() -> ViewRegistry:
    from sc_browser.views import (
        ClusterView,
        ExpressionView,
        FeatureCountView,
        DotplotView,
        HeatmapView,
        VolcanoPlotView,
        DatasetSummary,
    )

    registry = ViewRegistry()
    registry.register(ClusterView)
    registry.register(ExpressionView)
    registry.register(FeatureCountView)
    registry.register(DotplotView)
    registry.register(HeatmapView)
    registry.register(VolcanoPlotView)
    registry.register(DatasetSummary)
    return registry


def _choose_default_dataset(
    datasets: list[Dataset],
    global_config,
) -> Optional[Dataset]:
    """
    Pick the default dataset based on GlobalConfig.default_group if possible,
    otherwise fall back to the first dataset.

    If no datasets are available, returns None.
    """
    if not datasets:
        return None

    default_group = getattr(global_config, "default_group", None)
    if default_group:
        for ds in datasets:
            # Dataset should expose 'group' from its config; be defensive if missing
            if getattr(ds, "group", None) == default_group:
                return ds

    # Fallback: just use the first dataset
    return datasets[0]



def create_dash_app(config_root: Path | str = Path("config")) -> Dash:
    config_root = Path(config_root)

    # 1) config only (no datasets loaded)
    global_config, cfg_by_name = load_dataset_registry(config_root)
    if not cfg_by_name:
        raise RuntimeError("No dataset configs were loaded from config")

    # 2) lazy dataset maps
    dataset_cache: Dict[str, Dataset] = {}
    dataset_by_name = _LazyDatasetsByName(cfg_by_name=cfg_by_name, cache=dataset_cache)

    # dataset key mapping (uses config key if present, else name)
    name_by_key: Dict[str, str] = {}
    for name, cfg in cfg_by_name.items():
        key = getattr(cfg, "key", None) or name
        if key in name_by_key:
            raise RuntimeError(f"Duplicate dataset key '{key}' in config")
        name_by_key[key] = name

    dataset_by_key = _LazyDatasetsByKey(name_by_key=name_by_key, by_name=dataset_by_name)

    registry = _build_view_registry()

    # 3) choose + load ONLY default dataset
    default_dataset: Optional[Dataset] = None
    default_group = getattr(global_config, "default_group", None)

    if default_group:
        for cfg in cfg_by_name.values():
            if getattr(cfg, "group", None) == default_group:
                default_dataset = dataset_by_name[cfg.name]
                break

    if default_dataset is None:
        # deterministic fallback
        first_name = sorted(cfg_by_name.keys())[0]
        default_dataset = dataset_by_name[first_name]

    export_root = config_root / "exports"
    export_root.mkdir(exist_ok=True)

    export_service = ExportService(
        datasets_by_key=dataset_by_key,
        view_registry=registry,
        output_root=export_root,
    )

    # Keep ctx.datasets as "loaded datasets" if anything relies on it
    # (Only default loaded at boot)
    ctx = AppConfig(
        config_root=config_root,
        global_config=global_config,
        datasets=[default_dataset],
        dataset_by_name=dataset_by_name,
        dataset_by_key=dataset_by_key,
        default_dataset=default_dataset,
        registry=registry,
        export_service=export_service,
        enable_dataset_management=bool(os.getenv("ENABLE_DATASET_MANAGEMENT", "0") == "1"),
    )

    app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY])
    app.title = getattr(global_config, "ui_title", "Single-Cell Browser")

    app.layout = build_layout(ctx)
    register_explore_callbacks(app, ctx)
    register_dataset_import_callbacks(app, ctx)
    register_reports_callbacks(app, ctx)
    register_dataset_preview_callbacks(app, ctx)

    return app
