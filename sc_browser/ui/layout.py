from __future__ import annotations

from typing import List, TYPE_CHECKING, Optional

import dash_bootstrap_components as dbc
from dash import dcc, html

from sc_browser.core.dataset import Dataset
from sc_browser.core.view_registry import ViewRegistry

from .helpers import get_filter_dropdown_options

if TYPE_CHECKING:
    from .dash_app import AppConfig



def _build_navbar(datasets: List[Dataset], global_config, default_dataset: Dataset | None) -> dbc.Navbar:
    title = getattr(global_config, "ui_title", "sc-B++")
    subtitle = getattr(global_config, "subtitle", "Interactive Dataset Explorer")

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


def _build_view_and_label_panel(registry: ViewRegistry) -> dbc.Card:
    view_classes = registry.all_classes()
    view_options = [
        {"label": cls.label, "value": cls.id}
        for cls in view_classes
    ]
    default_view_id: str | None = view_classes[0].id if view_classes else None

    return dbc.Card(
        [
            dbc.CardHeader("View & label", className="fw-semibold"),
            dbc.CardBody(
                [
                    # ------------------------------
                    # 1) Saved / New view selector
                    # ------------------------------
                    html.Div(
                        [
                            html.Label("Saved view", className="form-label"),
                            dcc.Dropdown(
                                id="saved-figure-select",
                                options=[],  # populated from session-metadata
                                value="__new__",  # "New view" sentinel
                                placeholder="New view (start from current filters)",
                                clearable=False,
                                className="mb-1",
                            ),
                            html.Small(
                                "Pick 'New view' to work from your current filters, or choose a saved view "
                                "and click Load to restore it.",
                                className="text-muted",
                            ),
                        ],
                        className="mb-3",
                    ),

                    # ------------------------------
                    # 2) View type selector
                    # ------------------------------
                    html.Div(
                        [
                            html.Label("View type", className="form-label"),
                            dcc.Dropdown(
                                id="view-select",
                                options=view_options,
                                value=default_view_id,
                                clearable=False,
                                placeholder="Select view type",
                                className="mb-1",
                            ),
                            html.Small(
                                "Choose the kind of plot (clusters, expression, dotplot, etc.).",
                                className="text-muted",
                            ),
                        ],
                        className="mb-3",
                    ),

                    # ------------------------------
                    # 3) Figure label input
                    # ------------------------------
                    html.Div(
                        [
                            html.Label("Figure label", className="form-label"),
                            dcc.Input(
                                id="figure-label-input",
                                type="text",
                                placeholder="Figure label (optional)",
                                className="form-control",
                            ),
                        ],
                        className="mb-3",
                    ),

                    # ------------------------------
                    # 4) Load / Save buttons
                    # ------------------------------
                    html.Div(
                        [
                            dbc.Button(
                                "Load",
                                id="saved-figure-load-btn",
                                n_clicks=0,
                                color="secondary",
                                outline=True,
                                size="sm",
                                className="me-2",
                            ),
                            dbc.Button(
                                "Save",
                                id="save-figure-btn",
                                n_clicks=0,
                                color="primary",
                                size="sm",
                            ),
                        ],
                        className="mb-2",
                    ),

                    # ------------------------------
                    # 5) Status text for save operations
                    # ------------------------------
                    html.Div(id="save-figure-status", className="text-muted"),
                ]
            ),
        ],
        className="scb-viewlabel-card",
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
                    html.Div(
                        id="gene-filter-container",
                        children=[
                            html.Label("Gene(s)", className="form-label"),
                            dcc.Dropdown(
                                id="gene-select",
                                options=[],
                                multi=True,
                                placeholder="Type to search genes",
                                className="mb-3",
                            ),
                        ],
                    ),
                    html.Div(
                        id="embedding-filter-container",
                        children=[
                            html.Label("Embedding", className="form-label"),
                            dcc.Dropdown(
                                id="embedding-select",
                                options=emb_options,
                                value=(
                                    default_dataset.embedding_key
                                    if default_dataset.embedding_key in default_dataset.adata.obsm
                                    else None
                                ),
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
                    html.Div(
                        id="options-container",
                        children=[
                            dbc.Checklist(
                                id="options-checklist",
                                options=[
                                    {"label": " Split by condition", "value": "split_by_condition"},
                                    {"label": " 3D view", "value": "is_3d"},
                                ],
                                value=[],
                                switch=True,
                            ),
                            html.Hr(),
                        ],
                    ),
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


def _build_plot_panel() -> dbc.Card:
    return dbc.Card(
        [
            dbc.CardHeader(
                html.Div(
                    [
                        html.Strong("Plot"),
                    ],
                    className="d-flex align-items-center",
                ),
                className="p-2",
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
                                className="mt-2 ms-auto me-2",
                            ),
                            dcc.Download(id="download-data"),
                        ],
                        className="d-flex justify-content-end align-items-center",
                    ),
                ],
                className="scb-main-body",
            ),
        ],
        className="scb-maincard",
    )

def _build_dataset_import_panel() -> dbc.Container:
    """
    Datasets → Import & mapping page.

    - Upload/import .h5ad
    - Show current dataset status + summary
    - Map obs columns / embedding
    """
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
        ],
        className="scb-datasets-view",
    )


def _build_dataset_preview_panel() -> dbc.Container:
    """
    Dataset preview page:

    - Shows .obs preview (and later any other QC / summary)
    """
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
                    dbc.Col(obs_card, md=12, className="mt-3"),
                ],
                className="gx-3",
            ),
        ],
        className="scb-datasets-view",
    )


def _build_reports_panel() -> dbc.Container:
    """
    Reports page:

    - Shows how many figures are in the current session
    - Lists all figures (populated via callbacks in reports-figure-list)
    - Allows metadata_io of session metadata (+ images)
    - Allows import of a metadata JSON
    """
    header = dbc.Row(
        [
            dbc.Col(
                html.Div(
                    [
                        html.H4("Reports"),
                        html.Div(
                            id="reports-summary-text",
                            className="text-muted small",
                        ),
                    ]
                ),
                md=6,
            ),
            dbc.Col(
                html.Div(
                    [
                        dbc.Button(
                            "Export session",
                            id="reports-metadata_io-btn",
                            color="primary",
                            size="sm",
                            className="me-2",
                        ),
                        dcc.Download(id="reports-download-session"),

                        dcc.Upload(
                            id="reports-upload",
                            children=html.Div(
                                ["Import metadata JSON (drag & drop or ", html.A("select file"), ")"]
                            ),
                            multiple=False,
                            className="scb-upload border rounded p-2 text-center d-inline-block",
                        ),
                    ],
                    className="d-flex justify-content-end align-items-center",
                ),
                md=6,
            ),
        ],
        className="mt-3 mb-2",
    )

    banner = html.Div(
        id="reports-import-banner",
        className="small text-muted mb-2",
    )

    table_header = dbc.Row(
        [
            dbc.Col(html.Strong("Saved figures"), md=12),
        ],
        className="mt-3 mb-1",
    )

    # This div will be populated with a table/list of figures via callbacks
    figure_list = html.Div(
        id="reports-figure-list",
        className="scb-reports-figure-list",
    )

    return dbc.Container(
        fluid=True,
        children=[
            header,
            banner,
            html.Hr(),
            table_header,
            figure_list,
        ],
        className="scb-reports-view",
    )


def build_layout(ctx: "AppContext"):
    default_dataset = _choose_default_dataset(ctx.datasets, ctx.global_config)

    navbar = _build_navbar(ctx.datasets, ctx.global_config, default_dataset)

    status_bar = dbc.Alert(
        id="status-bar",
        color="light",
        className="my-2 py-2 scb-status-bar",
    )

    if default_dataset is None:
        filter_panel = dbc.Card(
            dbc.CardBody("No datasets configured. Please add a dataset in the Datasets tab."),
            className="scb-sidebar",
        )
        view_panel = dbc.Card(
            dbc.CardBody("No views available without a dataset."),
            className="scb-viewlabel-card",
        )
    else:
        filter_panel = _build_filter_panel(default_dataset)
        view_panel = _build_view_and_label_panel(ctx.registry)

    plot_panel = _build_plot_panel()
    dataset_import_panel = _build_dataset_import_panel()
    dataset_preview_panel = _build_dataset_preview_panel()
    reports_panel = _build_reports_panel()

    return dbc.Container(
        fluid=True,
        className="scb-root",
        children=[
            navbar,
            status_bar,

            # App-level stores
            dcc.Store(id="session-metadata", storage_type="session"),
            dcc.Store(id="active-session-id", storage_type="session"),
            dcc.Store(id="active-figure-id", storage_type="session"),
            dcc.Store(id="filter-state", storage_type="session"),
            dcc.Store(id="user-state", storage_type="local"),


            dcc.Tabs(
                id="page-tabs",
                value="explore",
                children=[
                    dcc.Tab(
                        label="Dataset Importer",
                        value="datasets",
                        children=[dataset_import_panel],
                    ),
                    dcc.Tab(
                        label="Dataset Preview",
                        value="dataset-preview",
                        children=[dataset_preview_panel],
                    ),
                    dcc.Tab(
                        label="Explore",
                        value="explore",
                        children=[
                            dbc.Row(
                                [
                                    dbc.Col(
                                        [
                                            view_panel,
                                            filter_panel,
                                        ],
                                        md=3,
                                        className="mt-3",
                                    ),
                                    dbc.Col(
                                        plot_panel,
                                        md=9,
                                        className="mt-3",
                                    ),
                                ],
                                className="gx-3",
                            ),
                        ],
                    ),
                    dcc.Tab(
                        label="Reports",
                        value="reports",
                        children=[reports_panel],
                    ),
                ],
                className="mt-2",
            ),
        ],
    )


def _choose_default_dataset(datasets: List[Dataset], global_config) -> Optional[Dataset]:
    if not datasets:
        return None

    default_group = getattr(global_config, "default_group", None)
    if default_group:
        for ds in datasets:
            if getattr(ds, "group", None) == default_group:
                return ds

    return datasets[0]
