from __future__ import annotations

import dash_bootstrap_components as dbc

from dash import dcc, html


def build_reports_panel() -> dbc.Container:
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