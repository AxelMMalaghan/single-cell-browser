from __future__ import annotations

import dash_bootstrap_components as dbc
from dash import dcc
from typing import List, Optional

from sc_browser.core.dataset import Dataset
from sc_browser.ui.layout.build_plot_panel import build_plot_panel
from sc_browser.ui.layout.build_dataset_import_panel import build_dataset_import_panel
from sc_browser.ui.layout.build_view_and_label_panel import build_view_and_label_panel
from sc_browser.ui.layout.build_navbar import build_navbar
from sc_browser.ui.layout.build_filter_panel import build_filter_panel
from sc_browser.ui.layout.build_reports_panel import build_reports_panel
from sc_browser.ui.layout.build_dataset_preview_panel import build_dataset_preview_panel


def build_layout(ctx: "AppContext"):
    default_dataset = _choose_default_dataset(ctx.datasets, ctx.global_config)

    navbar = build_navbar(ctx.datasets, ctx.global_config, default_dataset)

    if default_dataset is None:
        filter_panel = dbc.Card(
            dbc.CardBody("No datasets configured. Please add a dataset in the Datasets tab."),
            className="scb-sidebar",
        )
        view_panel = dbc.Card(
            dbc.CardBody("No views available without a dataset."),
            className="scb-viewlabel-card",
        )
    else:
        filter_panel = build_filter_panel(default_dataset)
        view_panel = build_view_and_label_panel(ctx.registry)

    plot_panel = build_plot_panel()
    dataset_import_panel = build_dataset_import_panel(ctx.datasets, default_dataset)
    dataset_preview_panel = build_dataset_preview_panel()
    reports_panel = build_reports_panel()

    return dbc.Container(
        fluid=True,
        className="scb-root",
        children=[
            navbar,

            # App-level stores
            dcc.Store(id="session-metadata", storage_type="session"),
            dcc.Store(id="active-session-id", storage_type="session"),
            dcc.Store(id="active-figure-id", storage_type="session"),
            dcc.Store(id="filter-state", storage_type="session"),
            dcc.Store(id="user-state", storage_type="local"),

            dcc.Tabs(
                id="page-tabs",
                value="explore",
                children=[
                    dcc.Tab(
                        label="Datasets",
                        value="datasets",
                        children=[dataset_import_panel, dataset_preview_panel],
                    ),
                    dcc.Tab(
                        label="Explore",
                        value="explore",
                        children=[
                            dbc.Row(
                                [
                                    dbc.Col(
                                        [
                                            view_panel,
                                            filter_panel,
                                        ],
                                        md=3,
                                        className="mt-3",
                                    ),
                                    dbc.Col(
                                        plot_panel,
                                        md=9,
                                        className="mt-3",
                                    ),
                                ],
                                className="gx-3",
                            ),
                        ],
                    ),
                    dcc.Tab(
                        label="Summary",
                        value="Summary",
                        children=[reports_panel],
                    ),
                ],
                className="mt-2",
            ),
        ],
    )


def _choose_default_dataset(datasets: List[Dataset], global_config) -> Optional[Dataset]:
    if not datasets:
        return None

    default_group = getattr(global_config, "default_group", None)
    if default_group:
        for ds in datasets:
            if getattr(ds, "group", None) == default_group:
                return ds

    return datasets[0]

