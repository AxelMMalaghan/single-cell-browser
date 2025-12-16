from __future__ import annotations

import dash_bootstrap_components as dbc

from dash import dcc, html

from sc_browser.ui.ids import IDs


def build_dataset_import_panel() -> dbc.Container:
    """
    Datasets → Import & mapping page.

    - Upload/import .h5ad
    - Show current dataset status + summary
    - Map obs columns / embedding
    """
    status_card = dbc.Card(
        [
            dbc.Toast(
                id=IDs.Control.DATASET_TOAST,  # create this ID
                header="Dataset",
                is_open=False,
                dismissable=True,
                duration=3500,
                icon="success",  # will be overridden by callback if needed
                style={
                    "position": "fixed",
                    "top": 80,
                    "right": 16,
                    "width": 360,
                    "zIndex": 2000,
                },
            ),
            dbc.CardHeader("Current dataset"),
            dbc.CardBody(
                [
                    # Changed from 'mb-1 fw-semibold' to 'form-label mb-1' to match right side labels
                    html.Div(id="dm-current-dataset", className="form-label mb-1"),

                    # Added 'small' class to match the font size of the labels (0.78rem)
                    html.Div(id="dm-status-text", className="mb-1 small"),
                    html.Div(id="dm-summary-text", className="text-muted mb-3 small"),

                    html.Hr(),

                    # Changed from html.H6 to html.Label with 'form-label' class
                    html.Label("Import dataset (.h5ad)", className="form-label mt-1"),

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

    html.Span(
        id=IDs.Control.DATASET_STATUS_BADGE,  # create this ID
        className="ms-2",
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