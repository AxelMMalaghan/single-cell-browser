from __future__ import annotations

import dash_bootstrap_components as dbc

from dash import dcc, html


def build_plot_panel() -> dbc.Card:
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