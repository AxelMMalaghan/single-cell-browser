from __future__ import annotations

from typing import List

import dash_bootstrap_components as dbc
from dash import html

from sc_browser.core.dataset import Dataset


def build_navbar(
    datasets: List[Dataset],
    global_config,
    default_dataset: Dataset | None,
) -> dbc.Navbar:
    navbar_image_src = getattr(
        global_config, "navbar_image_src", "/assets/hgtc_logo.png"
    )

    title = getattr(global_config, "ui_title", "Dashy")
    subtitle = getattr(global_config, "subtitle", "Interactive Dataset Explorer")

    return dbc.Navbar(
        dbc.Container(
            fluid=True,
            children=[
                # Left: logo + title
                html.Div(
                    className="d-flex align-items-center",
                    children=[
                        html.Img(
                            src=navbar_image_src,
                            alt="Logo",
                            style={"height": "70px"},
                            className="me-3",
                        ),
                        html.Div(
                            [
                                html.H2(title, className="mb-0"),
                                html.Small(
                                    subtitle,
                                    className="text-muted",
                                    id="navbar-subtitle",
                                ),
                            ],
                            className="d-flex flex-column justify-content-center",
                        ),
                    ],
                ),
            ],
        ),
        dark=False,
        className="shadow-sm scb-navbar",
    )
