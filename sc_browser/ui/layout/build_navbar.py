from __future__ import annotations

from typing import List

import dash_bootstrap_components as dbc
from dash import dcc, html

from sc_browser.core.dataset import Dataset


def build_navbar(
    datasets: List[Dataset],
    global_config,
    default_dataset: Dataset | None,
) -> dbc.Navbar:
    navbar_image_src = getattr(global_config, "navbar_image_src", "/assets/hgtc_logo.png")

    title = getattr(global_config, "ui_title", "Dashy")
    subtitle = getattr(global_config, "subtitle", "Interactive Dataset Explorer")

    if default_dataset is not None:
        default_name = default_dataset.name
    else:
        default_name = datasets[0].name if datasets else None

    dataset_options = [{"label": ds.name, "value": ds.name} for ds in datasets]

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
        ),
        dark=False,
        className="shadow-sm scb-navbar",
    )