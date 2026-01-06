from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING, Dict, Tuple

if TYPE_CHECKING:
    from sc_browser.views.cluster_view import ClusterView
    from sc_browser.views.dataset_summary_view import DatasetSummary
    from sc_browser.views.dot_plot_view import DotplotView
    from sc_browser.views.feature_count_view import FeatureCountView
    from sc_browser.views.gene_expression_view import ExpressionView
    from sc_browser.views.heat_map_view import HeatmapView
    from sc_browser.views.volcano_plot_view import VolcanoPlotView

_VIEW_IMPORTS: Dict[str, Tuple[str, str]] = {
    "ClusterView": ("sc_browser.views.cluster_view", "ClusterView"),
    "ExpressionView": ("sc_browser.views.gene_expression_view", "ExpressionView"),
    "FeatureCountView": ("sc_browser.views.feature_count_view", "FeatureCountView"),
    "DotplotView": ("sc_browser.views.dot_plot_view", "DotplotView"),
    "HeatmapView": ("sc_browser.views.heat_map_view", "HeatmapView"),
    "VolcanoPlotView": ("sc_browser.views.volcano_plot_view", "VolcanoPlotView"),
    "DatasetSummary": ("sc_browser.views.dataset_summary_view", "DatasetSummary"),
}

__all__ = [
    "ClusterView",
    "ExpressionView",
    "FeatureCountView",
    "DotplotView",
    "HeatmapView",
    "VolcanoPlotView",
    "DatasetSummary",
]


def __getattr__(name: str):
    target = _VIEW_IMPORTS.get(name)
    if target is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr = target
    module = import_module(module_name)
    return getattr(module, attr)
