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

    # -- 1. Action Buttons (Export / Import) --
    # We group these to put them in the top-right of the card body
    actions_div = html.Div(
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
                    ["Import JSON (", html.A("select"), ")"]
                ),
                multiple=False,
                className="scb-upload border rounded px-2 py-1 text-center d-inline-block small",
                style={"lineHeight": "1.5", "cursor": "pointer"}
            ),
        ],
        className="d-flex justify-content-end align-items-center",
    )

    # -- 2. Main Content inside a Card --
    main_card = dbc.Card(
        [
            dbc.CardHeader("Session Reports"),
            dbc.CardBody(
                [
                    # Top Row: Summary (left) + Actions (right)
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    # This div is populated by callback (e.g. "5 saved figures...")
                                    html.Div(
                                        id="reports-summary-text",
                                        className="text-muted small mb-1",
                                    ),
                                    # The import banner feedback
                                    html.Div(
                                        id="reports-import-banner",
                                        className="small text-success fw-bold",
                                    ),
                                ],
                                md=6,
                                className="d-flex flex-column justify-content-center"
                            ),
                            dbc.Col(actions_div, md=6),
                        ],
                        className="mb-3",
                    ),

                    html.Hr(),

                    # Sub-header using the standard 'form-label' style
                    html.Label("Saved figures", className="form-label"),

                    # The list of figures (table) goes here
                    html.Div(
                        id="reports-figure-list",
                        className="scb-reports-figure-list mt-2",
                    ),
                ]
            ),
        ],
        className="h-100 shadow-sm",
    )

    return dbc.Container(
        fluid=True,
        children=[
            dbc.Row(
                [
                    dbc.Col(main_card, md=12, className="mt-3"),
                ],
                className="gx-3",
            ),
        ],
        className="scb-reports-view",
    )