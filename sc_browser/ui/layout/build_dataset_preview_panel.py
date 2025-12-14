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