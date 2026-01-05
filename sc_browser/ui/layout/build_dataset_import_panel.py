from __future__ import annotations

from typing import List

import dash_bootstrap_components as dbc

from dash import dcc, html
from sc_browser.core.dataset import Dataset



def build_dataset_import_panel(
    datasets: List[Dataset],
    default_dataset: Dataset | None,
) -> dbc.Container:
    """
    Datasets → Import & mapping page.
    """

    dataset_options = [{"label": ds.name, "value": ds.name} for ds in datasets]

    if default_dataset is not None:
        default_name = default_dataset.name
    else:
        default_name = datasets[0].name if datasets else None

    status_card = dbc.Card(
        [

            dbc.CardHeader("Current dataset"),
            dbc.CardBody(
                [
                    html.Div(id="dm-current-dataset", className="form-label mb-1"),
                    html.Div(id="dm-status-text", className="mb-1 small"),
                    html.Div(id="dm-summary-text", className="text-muted mb-3 small"),

                    html.Hr(),

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

            html.Div(
                [
                    html.Div(
                        "Active Dataset",
                        className="navbar-dataset-title",
                    ),
                    html.Div(
                        "Selects the current dataset",
                        className="navbar-dataset-subtitle",
                    ),
                    dcc.Dropdown(
                        id="dataset-select",
                        options=dataset_options,
                        value=default_name,
                        clearable=False,
                        placeholder="Select dataset",
                        className="scb-dataset-dropdown mt-1",
                    ),

                ],
                className="ms-auto navbar-dataset-block",
                style={
                    "minWidth": "280px",
                    "maxWidth": "380px",
                    "marginRight": "24px",
                },
            )
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