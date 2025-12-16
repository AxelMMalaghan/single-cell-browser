from __future__ import annotations

import dash_bootstrap_components as dbc
from dash import dcc, html
from sc_browser.core.dataset import Dataset
from sc_browser.ui.helpers import get_filter_dropdown_options

def build_filter_panel(default_dataset: Dataset) -> dbc.Card:
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
                    # --- NEW SECTION: Dataset Info Header ---
                    html.Div([
                        html.H5(
                            default_dataset.name,
                            id="sidebar-dataset-name",
                            className="card-title"
                        ),
                        html.P(
                            f"{default_dataset.adata.n_obs} cells · {default_dataset.adata.n_vars} genes",
                            id="sidebar-dataset-meta",
                            className="card-subtitle text-muted mb-3"
                        ),
                        html.Hr(),
                    ]),
                    # ----------------------------------------

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
                        id="colour-scale-container",
                        children=[
                            html.Label("Colour scale", className="form-label"),
                            dcc.Dropdown(
                                id="colour-scale-select",
                                options=[
                                    {"label": "Viridis", "value": "viridis"},
                                    {"label": "Plasma", "value": "plasma"},
                                    {"label": "Reds", "value": "reds"},
                                    {"label": "Greys", "value": "Greys"},
                                    {"label": "Turbo", "value": "turbo"},
                                ],
                                value="viridis",
                                clearable=False,
                                className="mb-3",
                            ),
                        ],
                        className="mb-2",
                    ),
                    html.Div(
                        id="dim-filter-container",
                        children=[
                            html.Label("Dimension X", className="form-label"),
                            dcc.Dropdown(id="dim-x-select", className="mb-2"),

                            html.Label("Dimension Y", className="form-label"),
                            dcc.Dropdown(id="dim-y-select", className="mb-2"),

                            html.Div(
                                id="dim-z-container",
                                children=[
                                    html.Label("Dimension Z", className="form-label"),
                                    dcc.Dropdown(id="dim-z-select", className="mb-2"),
                                ]
                            ),
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
                ]
            ),
        ],
        className="scb-sidebar",
    )