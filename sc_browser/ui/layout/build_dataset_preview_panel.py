from __future__ import annotations

import dash_bootstrap_components as dbc

from dash import html


def build_dataset_preview_panel() -> dbc.Container:
    """
    Dataset preview page:

    - Shows .obs preview (and later any other QC / summary)
    """
    obs_card = dbc.Card(
        [
            # Use html.Strong inside a Div to match the 'Plot' panel header style exactly
            dbc.CardHeader(
                html.Div(
                    [
                        html.Strong("Observation Preview"),
                    ],
                    className="d-flex align-items-center",
                ),
                className="p-2",
            ),
            # Use 'scb-main-body' to match the padding of the main plot area
            dbc.CardBody(
                html.Div(id="dm-obs-preview"),
                className="scb-main-body",
            ),
        ],
        # 'scb-maincard' provides the border, shadow, and background
        className="scb-maincard mt-3",
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
