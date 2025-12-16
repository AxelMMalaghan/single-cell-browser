from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Dict, Optional

import dash_bootstrap_components as dbc
from dash import Dash

from sc_browser.config.dataset_loader import load_dataset_registry
from sc_browser.core.dataset import Dataset
from sc_browser.core.view_registry import ViewRegistry
from sc_browser.services.dataset_service import DatasetKeyManager, DatasetManager
from sc_browser.services.metadata_export_service import ExportService
from sc_browser.services.storage import LocalFileSystemStorage
from sc_browser.ui.callbacks.callbacks_dataset_import import (
    register_dataset_import_callbacks,
)
from sc_browser.ui.callbacks.callbacks_dataset_preview import (
    register_dataset_preview_callbacks,
)
from sc_browser.ui.callbacks.callbacks_filters import register_filter_callbacks
from sc_browser.ui.callbacks.callbacks_io import register_io_callbacks
from sc_browser.ui.callbacks.callbacks_render import register_render_callbacks
from sc_browser.ui.callbacks.callbacks_reports import register_reports_callbacks
from sc_browser.ui.callbacks.callbacks_sync import register_sync_callbacks
from sc_browser.ui.config import AppConfig
from sc_browser.ui.layout.build_layout import build_layout

logger = logging.getLogger(__name__)


def _build_view_registry() -> ViewRegistry:
    from sc_browser.views import (
        ClusterView,
        DatasetSummary,
        DotplotView,
        ExpressionView,
        FeatureCountView,
        HeatmapView,
        VolcanoPlotView,
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


def create_dash_app(config_root: Path | str = Path("config")) -> Dash:
    config_root = Path(config_root)

    # 1) Load Config
    global_config, cfg_by_name = load_dataset_registry(config_root)
    if not cfg_by_name:
        raise RuntimeError("No dataset configs were loaded from config")

    # 2) Initialize Service Layer
    dataset_manager = DatasetManager(cfg_by_name)

    # Build key mapping (uses config key if present, else name)
    name_by_key: Dict[str, str] = {}
    for name, cfg in cfg_by_name.items():
        key = getattr(cfg, "key", None) or name
        if key in name_by_key:
            raise RuntimeError(f"Duplicate dataset key '{key}' in config")
        name_by_key[key] = name

    dataset_key_manager = DatasetKeyManager(dataset_manager, name_by_key)
    registry = _build_view_registry()

    # 3) Choose Default Dataset
    default_dataset: Optional[Dataset] = None
    default_group = getattr(global_config, "default_group", None)

    if default_group:
        for cfg in cfg_by_name.values():
            if getattr(cfg, "group", None) == default_group:
                # Trigger lazy load for the default
                default_dataset = dataset_manager.get(cfg.name)
                break

    if default_dataset is None:
        # Fallback: first available
        first_name = sorted(cfg_by_name.keys())[0]
        default_dataset = dataset_manager.get(first_name)

    # 4) Export Service & Storage
    export_root = config_root / "exports"

    # Initialize storage backend (local fs)
    storage_backend = LocalFileSystemStorage(export_root)

    export_service = ExportService(
        datasets_by_key=dataset_key_manager,
        view_registry=registry,
        storage=storage_backend,
    )

    # 5) App Context
    ctx = AppConfig(
        config_root=config_root,
        global_config=global_config,
        datasets=[default_dataset],  # Keep compatible with list-based access if needed
        dataset_by_name=dataset_manager,
        dataset_by_key=dataset_key_manager,
        default_dataset=default_dataset,
        registry=registry,
        export_service=export_service,
        enable_dataset_management=bool(os.getenv("ENABLE_DATASET_MANAGEMENT", "0") == "1"),
    )

    app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY])
    app.title = getattr(global_config, "ui_title", "Single-Cell Browser")

    app.layout = build_layout(ctx)

    # Register callbacks
    register_sync_callbacks(app, ctx)
    register_io_callbacks(app, ctx)
    register_render_callbacks(app, ctx)
    register_filter_callbacks(app, ctx)
    register_dataset_import_callbacks(app, ctx)
    register_reports_callbacks(app, ctx)
    register_dataset_preview_callbacks(app, ctx)

    return app