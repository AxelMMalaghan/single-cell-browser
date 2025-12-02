from __future__ import annotations

from typing import Dict, List, Tuple

from dash import Dash, dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
import pandas as pd

from sc_browser.config.io import load_datasets
from sc_browser.core.state import FilterState
from sc_browser.core.view_registry import ViewRegistry
from sc_browser.views import (
    ClusterView,
    ExpressionView,
    FeatureCountView,
    Dotplot,
    HeatmapView,
    VolcanoPlotView, DatasetSummary,
)

def _build_view_registry() -> ViewRegistry:
    registry = ViewRegistry()
    registry.register(ClusterView)
    registry.register(ExpressionView)
    registry.register(FeatureCountView)
    registry.register(Dotplot)
    registry.register(HeatmapView)
    registry.register(VolcanoPlotView)
    registry.register(DatasetSummary)
    return registry


def _build_navbar(datasets) -> dbc.Navbar:
    """
    Top navbar: title on the left, dataset selector on the right.
    """
    return dbc.Navbar(
        dbc.Container(
            fluid=True,
            children=[
                # Title + subtitle on the left
                html.Div(
                    [
                        html.H2("Single-Cell Browser", className="mb-0"),
                        html.Small(
                            "Interactive Dataset Explorer",
                            className="text-muted",
                            id="navbar-subtitle",
                        ),
                    ],
                    className="d-flex flex-column justify-content-center",
                ),
                # Dataset selector on the right
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
        color="light",  # theme primary colour
        dark=False,        # white text on coloured navbar
        className="shadow-sm scb-navbar",
    )


def _build_filter_panel(default_dataset) -> dbc.Card:
    """
    Left sidebar with filters (cluster, condition, genes, options)
    """
    return dbc.Card(
        [
            dbc.CardHeader("Filters", className="fw-semibold"),
            dbc.CardBody(
                [
                    html.Label("Filter clusters", className="form-label"),
                    dcc.Dropdown(
                        id="cluster-select",
                        options=[
                            {"label": c, "value": c}
                            for c in sorted(default_dataset.clusters.unique())
                        ],
                        multi=True,
                        placeholder="All clusters",
                        className="mb-3",
                    ),
                    html.Label("Filter conditions", className="form-label"),
                    dcc.Dropdown(
                        id="condition-select",
                        options=[
                            {"label": c, "value": c}
                            for c in sorted(default_dataset.conditions.unique())
                        ],
                        multi=True,
                        placeholder="All conditions",
                        className="mb-3",
                    ),
                    html.Label("Gene(s)", className="form-label"),
                    dcc.Dropdown(
                        id="gene-select",
                        options=[
                            {"label": g, "value": g}
                            for g in sorted(default_dataset.genes)
                        ],
                        multi=True,
                        placeholder="Select gene(s)",
                        className="mb-3",
                    ),
                    html.Hr(),
                    dbc.Checklist(
                        id="options-checklist",
                        options=[
                            {
                                "label": " Split by condition",
                                "value": "split_by_condition",
                            },
                        ],
                        value=[],
                        switch=True,
                    ),
                    html.Hr(),
                    # Tiny dataset summary
                    html.Div(
                        [
                            html.Div(
                                f"{default_dataset.name}",
                                className="fw-semibold",
                            ),
                            html.Div(
                                [
                                    html.Span(
                                        f"{default_dataset.adata.n_obs} cells · "
                                        f"{default_dataset.adata.n_vars} genes",
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
                        dcc.Tab(label=cls.label, value=cls.id)
                        for cls in registry.all_classes()
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


# ----------------------------------------------------------------------
# App factory
# ----------------------------------------------------------------------


def create_dash_app() -> Dash:
    """
    Build and wire a complete Dash app around the core domain objects
    (Dataset, FilterState, ViewRegistry, Views).
    """
    # ----- Load datasets -----
    global_config, datasets = load_datasets("config/datasets.json")
    if not datasets:
        raise RuntimeError("No datasets were loaded from config/datasets.json")

    dataset_by_name: Dict[str, object] = {ds.name: ds for ds in datasets}
    default_dataset = datasets[0]

    # ----- View registry -----
    registry = _build_view_registry()

    # ----- Create Dash app with a coloured theme -----
    app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY])

    navbar = _build_navbar(datasets)
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

    # ------------------------------------------------------------------
    # Callbacks
    # ------------------------------------------------------------------

    @app.callback(
        Output("cluster-select", "options"),
        Output("condition-select", "options"),
        Output("gene-select", "options"),
        Input("dataset-select", "value"),
    )
    def update_filters(dataset_name: str):
        ds = dataset_by_name[dataset_name]
        cluster_options = [{"label": c, "value": c} for c in sorted(ds.clusters.unique())]
        condition_options = [{"label": c, "value": c} for c in sorted(ds.conditions.unique())]
        gene_options = [{"label": g, "value": g} for g in sorted(ds.genes)]
        return cluster_options, condition_options, gene_options

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
        fig = view.render_figure(data, state)
        return fig


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

        # If view returns non-tabular data, just bail for now
        if not isinstance(data, pd.DataFrame) or data.empty:
            return None

        filename = f"{view_id}_{dataset_name.replace(' ', '_')}.csv"

        # dcc.send_data_frame will call df.to_csv(...) under the hood
        return dcc.send_data_frame(data.to_csv, filename, index=False)

    return app