from __future__ import annotations

from pathlib import Path
from typing import Dict

import dash_bootstrap_components as dbc
from dash import Dash

from sc_browser.config.io import load_datasets
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


def create_dash_app(config_root: Path | str = Path("config")) -> Dash:
    config_root = Path(config_root)

    global_config, datasets = load_datasets(config_root)
    if not datasets:
        raise RuntimeError("No datasets were loaded from config")

    dataset_by_name: Dict[str, Dataset] = {ds.name: ds for ds in datasets}
    registry = _build_view_registry()

    ctx = AppContext(
        config_root=config_root,
        global_config=global_config,
        datasets=datasets,
        dataset_by_name=dataset_by_name,
        registry=registry,
    )

    app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY])

    app.layout = build_layout(ctx)

    register_explore_callbacks(app, ctx)
    register_dataset_callbacks(app, ctx)

    return app