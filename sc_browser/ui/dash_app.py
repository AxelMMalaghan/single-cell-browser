from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

from dash import Dash, dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
import pandas as pd

from sc_browser.config import GlobalConfig
from sc_browser.config.io import load_datasets
from sc_browser.core.state import FilterState
from sc_browser.core.view_registry import ViewRegistry
from sc_browser.views import (
    ClusterView,
    ExpressionView,
    FeatureCountView,
    DotplotView,
    HeatmapView,
    VolcanoPlotView,
    DatasetSummary,
)

def get_filter_dropdown_options(dataset) -> Tuple[List[dict], List[dict], List[dict]]:
    cluster_options = [{"label": c, "value": c} for c in sorted(dataset.clusters.unique())]
    condition_options = [{"label": c, "value": c} for c in sorted(dataset.conditions.unique())]
    gene_options = [{"label": g, "value": g} for g in sorted(dataset.genes)]
    return cluster_options, condition_options, gene_options

def _build_view_registry() -> ViewRegistry:
    registry = ViewRegistry()
    registry.register(ClusterView)
    registry.register(ExpressionView)
    registry.register(FeatureCountView)
    registry.register(DotplotView)
    registry.register(HeatmapView)
    registry.register(VolcanoPlotView)
    registry.register(DatasetSummary)
    return registry

def _build_navbar(datasets, global_config) -> dbc.Navbar:
    title = getattr(global_config, "ui_title", "Single-Cell Browser")
    subtitle = getattr(global_config, "subtitle", "Interactive Dataset Explorer")
    return dbc.Navbar(
        dbc.Container(
            fluid=True,
            children=[
                html.Div(
                    [
                        html.H2(title, className="mb-0"),
                        html.Small(subtitle, className="text-muted", id="navbar-subtitle"),
                    ],
                    className="d-flex flex-column justify-content-center",
                ),
                html.Div(
                    dcc.Dropdown(
                        id="dataset-select",
                        options=[{"label": ds.name, "value": ds.name} for ds in datasets],
                        value=datasets[0].name if datasets else None,
                        clearable=False,
                        placeholder="Select dataset",
                        style={"minWidth": "260px"},
                        className="scb-dataset-dropdown",
                    ),
                    className="ms-auto d-flex align-items-center",
                    style={"width": "300px"},
                ),
            ],
        ),
        color="light",
        dark=False,
        className="shadow-sm scb-navbar",
    )

def _build_filter_panel(default_dataset) -> dbc.Card:
    cluster_options, condition_options, gene_options = get_filter_dropdown_options(default_dataset)
    return dbc.Card(
        [
            dbc.CardHeader("Filters", className="fw-semibold"),
            dbc.CardBody(
                [
                    html.Label("Filter clusters", className="form-label"),
                    dcc.Dropdown(
                        id="cluster-select",
                        options=cluster_options,
                        multi=True,
                        placeholder="All clusters",
                        className="mb-3",
                    ),
                    html.Label("Filter conditions", className="form-label"),
                    dcc.Dropdown(
                        id="condition-select",
                        options=condition_options,
                        multi=True,
                        placeholder="All conditions",
                        className="mb-3",
                    ),
                    html.Label("Gene(s)", className="form-label"),
                    dcc.Dropdown(
                        id="gene-select",
                        options=gene_options,
                        multi=True,
                        placeholder="Select gene(s)",
                        className="mb-3",
                    ),
                    html.Hr(),
                    dbc.Checklist(
                        id="options-checklist",
                        options=[{"label": " Split by condition", "value": "split_by_condition"}],
                        value=[],
                        switch=True,
                    ),
                    html.Hr(),
                    html.Div(
                        [
                            html.Div(default_dataset.name, className="fw-semibold"),
                            html.Div(
                                [
                                    html.Span(
                                        f"{default_dataset.adata.n_obs} cells · {default_dataset.adata.n_vars} genes",
                                        className="text-muted",
                                    ),
                                ],
                                className="small",
                            ),
                        ],
                        className="scb-dataset-summary mt-1",
                    ),
                ]
            ),
        ],
        className="scb-sidebar",
    )

def _build_plot_panel(registry: ViewRegistry) -> dbc.Card:
    return dbc.Card(
        [
            dbc.CardHeader(
                dcc.Tabs(
                    id="view-tabs",
                    value=ClusterView.id,
                    children=[
                        dcc.Tab(label=cls.label, value=cls.id) for cls in registry.all_classes()
                    ],
                    className="scb-tabs",
                ),
                className="p-0",
            ),
            dbc.CardBody(
                [
                    dcc.Loading(
                        id="main-graph-loading",
                        type="default",
                        children=dcc.Graph(
                            id="main-graph",
                            style={"height": "650px"},
                            config={"responsive": True},
                        ),
                    ),
                    html.Div(
                        [
                            dbc.Button(
                                "Download data (CSV)",
                                id="download-data-btn",
                                color="secondary",
                                size="sm",
                                className="mt-2 me-2",
                            ),
                            dcc.Download(id="download-data"),
                        ],
                        className="d-flex justify-content-end",
                    ),
                ],
                className="scb-main-body",
            ),
        ],
        className="scb-maincard",
    )

def create_dash_app() -> Dash:

    config_root = Path("config")

    global_config, datasets = load_datasets(config_root)

    if not datasets:
        raise RuntimeError("No datasets were loaded from config/demo1dataset.json")

    dataset_by_name: Dict[str, object] = {ds.name: ds for ds in datasets}
    default_dataset = datasets[0]

    registry = _build_view_registry()

    app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY])

    navbar = _build_navbar(datasets, global_config)
    filter_panel = _build_filter_panel(default_dataset)
    plot_panel = _build_plot_panel(registry)

    app.layout = dbc.Container(
        fluid=True,
        className="scb-root",
        children=[
            navbar,
            dbc.Row(
                [
                    dbc.Col(filter_panel, md=3, className="mt-3"),
                    dbc.Col(plot_panel, md=9, className="mt-3"),
                ],
                className="gx-3",
            ),
        ],
    )

    @app.callback(
        Output("cluster-select", "options"),
        Output("condition-select", "options"),
        Output("gene-select", "options"),
        Input("dataset-select", "value"),
    )
    def update_filters(dataset_name: str):
        return get_filter_dropdown_options(dataset_by_name[dataset_name])

    @app.callback(
        Output("main-graph", "figure"),
        Input("view-tabs", "value"),
        Input("dataset-select", "value"),
        Input("cluster-select", "value"),
        Input("condition-select", "value"),
        Input("gene-select", "value"),
        Input("options-checklist", "value"),
    )
    def update_main_graph(
        view_id: str,
        dataset_name: str,
        clusters: List[str] | None,
        conditions: List[str] | None,
        genes: List[str] | None,
        options: List[str] | None,
    ):
        ds = dataset_by_name[dataset_name]
        state = FilterState(
            genes=genes or [],
            clusters=clusters or [],
            conditions=conditions or [],
            split_by_condition="split_by_condition" in (options or []),
        )
        view = registry.create(view_id, ds)
        data = view.compute_data(state)
        return view.render_figure(data, state)

    @app.callback(
        Output("download-data", "data"),
        Input("download-data-btn", "n_clicks"),
        State("view-tabs", "value"),
        State("dataset-select", "value"),
        State("cluster-select", "value"),
        State("condition-select", "value"),
        State("gene-select", "value"),
        State("options-checklist", "value"),
        prevent_initial_call=True,
    )
    def download_current_data(
        n_clicks,
        view_id: str,
        dataset_name: str,
        clusters,
        conditions,
        genes,
        options,
    ):
        ds = dataset_by_name[dataset_name]
        state = FilterState(
            genes=genes or [],
            clusters=clusters or [],
            conditions=conditions or [],
            split_by_condition="split_by_condition" in (options or []),
        )
        view = registry.create(view_id, ds)
        data = view.compute_data(state)
        if not isinstance(data, pd.DataFrame) or data.empty:
            return None
        filename = f"{view_id}_{dataset_name.replace(' ', '_')}.csv"
        return dcc.send_data_frame(data.to_csv, filename, index=False)

    return app
