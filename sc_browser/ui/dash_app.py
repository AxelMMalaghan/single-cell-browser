from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import dash_bootstrap_components as dbc
from dash import Dash

from sc_browser.config.loader import load_datasets
from sc_browser.core.view_registry import ViewRegistry
from sc_browser.core.dataset import Dataset

from .context import AppContext
from .layout import build_layout
from .callbacks_explore import register_explore_callbacks
from .callbacks_datasets import register_dataset_callbacks


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

    # Load config + datasets using the new IO layer
    global_config, datasets = load_datasets(config_root)
    if not datasets:
        raise RuntimeError("No datasets were loaded from config")

    dataset_by_name: Dict[str, Dataset] = {ds.name: ds for ds in datasets}
    registry = _build_view_registry()

    # Decide default dataset once, based on global_config
    default_dataset = _choose_default_dataset(datasets, global_config)

    ctx = AppContext(
        config_root=config_root,
        global_config=global_config,
        datasets=datasets,
        dataset_by_name=dataset_by_name,
        default_dataset=default_dataset,
        registry=registry,
    )

    app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY])

    # Drive browser window title from config
    app.title = getattr(global_config, "ui_title", "Single-Cell Browser")

    # Layout + callbacks both receive the same context
    app.layout = build_layout(ctx)

    register_explore_callbacks(app, ctx)
    register_dataset_callbacks(app, ctx)

    return app