from __future__ import annotations

from importlib import import_module
from typing import Dict, Tuple

_VIEW_IMPORTS: Dict[str, Tuple[str, str]] = {
    "ClusterView": ("sc_browser.views.cluster_view", "ClusterView"),
    "ExpressionView": ("sc_browser.views.gene_expression_view", "ExpressionView"),
    "FeatureCountView": ("sc_browser.views.feature_count_view", "FeatureCountView"),
    "DotplotView": ("sc_browser.views.dot_plot_view", "DotplotView"),
    "HeatmapView": ("sc_browser.views.heat_map_view", "HeatmapView"),
    "VolcanoPlotView": ("sc_browser.views.volcano_plot_view", "VolcanoPlotView"),
    "DatasetSummary": ("sc_browser.views.dataset_summary_view", "DatasetSummary"),
}

__all__ = list(_VIEW_IMPORTS.keys())


def __getattr__(name: str):
    target = _VIEW_IMPORTS.get(name)
    if target is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr = target
    module = import_module(module_name)
    return getattr(module, attr)
