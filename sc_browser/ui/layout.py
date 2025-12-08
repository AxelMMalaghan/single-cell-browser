# sc_browser/ui/layout.py

from __future__ import annotations

from typing import List, TYPE_CHECKING, Optional

import dash_bootstrap_components as dbc
from dash import dcc, html

from sc_browser.core.dataset import Dataset
from sc_browser.core.view_registry import ViewRegistry

from .helpers import get_filter_dropdown_options

if TYPE_CHECKING:
    from .dash_app import AppContext


def _build_navbar(datasets: List[Dataset], global_config, default_dataset: Dataset | None) -> dbc.Navbar:
    title = getattr(global_config, "ui_title", "scB++")
    subtitle = getattr(global_config, "subtitle", "Interactive Dataset Explorer")

    # Decide which dataset is initially selected in the dropdown
    if default_dataset is not None:
        default_name = default_dataset.name
    else:
        default_name = datasets[0].name if datasets else None

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
                        value=default_name,
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

def _build_filter_panel(default_dataset: Dataset) -> dbc.Card:
    (
        cluster_options,
        condition_options,
        sample_options,
        celltype_options,
        emb_options,
    ) = get_filter_dropdown_options(default_dataset)

    return dbc.Card(
        [
            dbc.CardHeader("Filters", className="fw-semibold"),
            dbc.CardBody(
                [
                    # Cluster filter
                    html.Div(
                        id="cluster-filter-container",
                        children=[
                            html.Label("Filter clusters", className="form-label"),
                            dcc.Dropdown(
                                id="cluster-select",
                                options=cluster_options,
                                multi=True,
                                placeholder="All clusters",
                                className="mb-3",
                            ),
                        ],
                    ),
                    # Condition filter
                    html.Div(
                        id="condition-filter-container",
                        children=[
                            html.Label("Filter conditions", className="form-label"),
                            dcc.Dropdown(
                                id="condition-select",
                                options=condition_options,
                                multi=True,
                                placeholder="All conditions",
                                className="mb-3",
                            ),
                        ],
                    ),
                    # Sample filter
                    html.Div(
                        id="sample-filter-container",
                        children=[
                            html.Label("Filter samples", className="form-label"),
                            dcc.Dropdown(
                                id="sample-select",
                                options=sample_options,
                                multi=True,
                                placeholder="All samples",
                                className="mb-3",
                            ),
                        ],
                    ),
                    # Cell type filter
                    html.Div(
                        id="celltype-filter-container",
                        children=[
                            html.Label("Filter cell types", className="form-label"),
                            dcc.Dropdown(
                                id="celltype-select",
                                options=celltype_options,
                                multi=True,
                                placeholder="All cell types",
                                className="mb-3",
                            ),
                        ],
                    ),
                    # Gene selector
                    html.Div(
                        id="gene-filter-container",
                        children=[
                            html.Label("Gene(s)", className="form-label"),
                            dcc.Dropdown(
                                id="gene-select",
                                options=[],  # populated via server-side search callback
                                multi=True,
                                placeholder="Type to search genes",
                                className="mb-3",
                            ),
                        ],
                    ),
                    # Embedding selector (e.g. PCA / TSNE / UMAP)
                    html.Div(
                        id="embedding-filter-container",
                        children=[
                            html.Label("Embedding", className="form-label"),
                            dcc.Dropdown(
                                id="embedding-select",
                                options=emb_options,
                                value=default_dataset.embedding_key
                                if default_dataset.embedding_key in default_dataset.adata.obsm
                                else None,
                                clearable=False,
                                placeholder="Select embedding",
                                className="mb-3",
                            ),
                        ],
                    ),
                    html.Div(
                        id="dim-filter-container",
                        children=[
                            html.Label("Dimension X", className="form-label"),
                            dcc.Dropdown(id="dim-x-select", className="mb-2"),
                            html.Label("Dimension Y", className="form-label"),
                            dcc.Dropdown(id="dim-y-select", className="mb-2"),
                            html.Label("Dimension Z", className="form-label"),
                            dcc.Dropdown(id="dim-z-select", className="mb-2"),
                        ],
                        className="mb-3",
                    ),

                    html.Hr(),
                    dbc.Checklist(
                        id="options-checklist",
                        options=[{"label": " Split by condition", "value": "split_by_condition"},
                                 {"label": "3D view", "value": "is_3d"}],
                        value=[],
                        switch=True,
                    ),
                    html.Hr(),
                    html.Div(
                        [
                            html.Div(
                                default_dataset.name,
                                id="sidebar-dataset-name",
                                className="fw-semibold",
                            ),
                            html.Div(
                                [
                                    html.Span(
                                        f"{default_dataset.adata.n_obs} cells · {default_dataset.adata.n_vars} genes",
                                        id="sidebar-dataset-meta",
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
                    value="cluster",  # will be overridden by first registered view id in callbacks if needed
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


def _build_dataset_manager_panel() -> dbc.Container:
    status_card = dbc.Card(
        [
            dbc.CardHeader("Current dataset"),
            dbc.CardBody(
                [
                    html.Div(id="dm-current-dataset", className="mb-1 fw-semibold"),
                    html.Div(id="dm-status-text", className="mb-1"),
                    html.Div(id="dm-summary-text", className="text-muted mb-3"),
                    html.Hr(),
                    html.H6("Import dataset (.h5ad)", className="mt-1"),
                    dcc.Upload(
                        id="dm-upload",
                        children=html.Div(
                            ["Drag and drop or ", html.A("select a .h5ad file")]
                        ),
                        multiple=False,
                        className="scb-upload border rounded p-2 text-center",
                    ),
                    html.Div(id="dm-import-status", className="small text-muted mt-2"),
                ]
            ),
        ],
        className="h-100",
    )

    mapping_card = dbc.Card(
        [
            dbc.CardHeader("Mapping"),
            dbc.CardBody(
                [
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    html.Label("Cluster column", className="form-label"),
                                    dcc.Dropdown(id="dm-cluster-key", className="mb-2"),
                                ],
                                md=6,
                            ),
                            dbc.Col(
                                [
                                    html.Label("Condition column", className="form-label"),
                                    dcc.Dropdown(id="dm-condition-key", className="mb-2"),
                                ],
                                md=6,
                            ),
                        ],
                        className="gx-3",
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    html.Label("Sample column", className="form-label"),
                                    dcc.Dropdown(id="dm-sample-key", className="mb-2"),
                                ],
                                md=6,
                            ),
                            dbc.Col(
                                [
                                    html.Label("Cell type column", className="form-label"),
                                    dcc.Dropdown(id="dm-celltype-key", className="mb-2"),
                                ],
                                md=6,
                            ),
                        ],
                        className="gx-3",
                    ),
                    html.Label("Embedding key", className="form-label mt-2"),
                    dcc.Dropdown(id="dm-embedding-key", className="mb-3"),
                    dbc.Button(
                        "Save mapping",
                        id="dm-save-btn",
                        color="primary",
                        size="sm",
                        className="mt-1",
                    ),
                    html.Span(id="dm-save-status", className="ms-2 small text-muted"),
                ]
            ),
        ],
        className="h-100",
    )

    obs_card = dbc.Card(
        [
            dbc.CardHeader("Obs preview"),
            dbc.CardBody(html.Div(id="dm-obs-preview"), className="p-2"),
        ],
        className="mt-3",
    )

    return dbc.Container(
        fluid=True,
        children=[
            dbc.Row(
                [
                    dbc.Col(status_card, md=4, className="mt-3"),
                    dbc.Col(mapping_card, md=8, className="mt-3"),
                ],
                className="gx-3",
            ),
            dbc.Row([dbc.Col(obs_card, md=12)], className="gx-3"),
        ],
        className="scb-datasets-view",
    )


def build_layout(ctx: "AppContext"):
    # Decide which dataset is the default for the whole UI
    default_dataset = _choose_default_dataset(ctx.datasets, ctx.global_config)

    navbar = _build_navbar(ctx.datasets, ctx.global_config, default_dataset)

    if default_dataset is None:
        # Graceful empty state if no datasets are configured
        filter_panel = dbc.Card(
            dbc.CardBody("No datasets configured. Please add a dataset in the Datasets tab."),
            className="scb-sidebar",
        )
    else:
        filter_panel = _build_filter_panel(default_dataset)

    plot_panel = _build_plot_panel(ctx.registry)
    dataset_manager = _build_dataset_manager_panel()

    return dbc.Container(
        fluid=True,
        className="scb-root",
        children=[
            navbar,
            dcc.Tabs(
                id="page-tabs",
                value="explore",
                children=[
                    dcc.Tab(
                        label="Explore",
                        value="explore",
                        children=[
                            dbc.Row(
                                [
                                    dbc.Col(filter_panel, md=3, className="mt-3"),
                                    dbc.Col(plot_panel, md=9, className="mt-3"),
                                ],
                                className="gx-3",
                            ),
                        ],
                    ),
                    dcc.Tab(
                        label="Datasets",
                        value="datasets",
                        children=[dataset_manager],
                    ),
                ],
                className="mt-2",
            ),
        ],
    )


def _choose_default_dataset(datasets: List[Dataset], global_config) -> Optional[Dataset]:
    """
    Pick the default dataset based on global_config.default_group if possible,
    otherwise fall back to the first dataset.

    If no datasets are available, returns None.
    """
    if not datasets:
        return None

    default_group = getattr(global_config, "default_group", None)
    if default_group:
        for ds in datasets:
            # Dataset should expose 'group' via config; be defensive if missing
            if getattr(ds, "group", None) == default_group:
                return ds

    # Fallback: just use the first dataset
    return datasets[0]
