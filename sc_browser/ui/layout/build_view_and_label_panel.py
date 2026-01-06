from __future__ import annotations

import dash_bootstrap_components as dbc
from dash import dcc, html

from sc_browser.core.view_registry import ViewRegistry


def build_view_and_label_panel(registry: ViewRegistry) -> dbc.Card:
    view_classes = registry.all_classes()
    view_options = [{"label": cls.label, "value": cls.id} for cls in view_classes]
    default_view_id: str | None = view_classes[0].id if view_classes else None

    return dbc.Card(
        [
            dbc.CardHeader("View & label", className="fw-semibold"),
            dbc.CardBody(
                [
                    # ------------------------------
                    # 1) Saved / New view selector
                    # ------------------------------
                    html.Div(
                        [
                            html.Label("Selected view", className="form-label"),
                            dbc.Select(
                                id="saved-figure-select",
                                options=[],  # populated from session-metadata
                                value="__new__",  # "New view" sentinel
                                className="scb-saved-view-select mb-1",
                            ),
                            html.Small(
                                "Pick 'New view' to work from your current filters, or choose a saved view "
                                "and click Load to restore it.",
                                className="text-muted",
                            ),
                        ],
                        className="mb-3",
                    ),
                    # ------------------------------
                    # 2) View type selector
                    # ------------------------------
                    html.Div(
                        [
                            html.Label("View type", className="form-label"),
                            dcc.Dropdown(
                                id="view-select",
                                options=view_options,
                                value=default_view_id,
                                clearable=False,
                                placeholder="Select view type",
                                className="mb-1",
                            ),
                            html.Small(
                                "Choose the kind of plot (clusters, expression, dotplot, etc.).",
                                className="text-muted",
                            ),
                        ],
                        className="mb-3",
                    ),
                    # ------------------------------
                    # 3) Figure label input
                    # ------------------------------
                    html.Div(
                        [
                            html.Label("Figure label", className="form-label"),
                            dcc.Input(
                                id="figure-label-input",
                                type="text",
                                placeholder="Figure label (optional)",
                                className="form-control",
                            ),
                        ],
                        className="mb-3",
                    ),
                    # ------------------------------
                    # 4) Load / Save buttons
                    # ------------------------------
                    html.Div(
                        [
                            dbc.Button(
                                "Load",
                                id="saved-figure-load-btn",
                                n_clicks=0,
                                color="secondary",
                                outline=True,
                                size="sm",
                                className="me-2",
                            ),
                            dbc.Button(
                                "Save",
                                id="save-figure-btn",
                                n_clicks=0,
                                color="primary",
                                size="sm",
                            ),
                        ],
                        className="mb-2",
                    ),
                    # ------------------------------
                    # 5) Status text for save operations
                    # ------------------------------
                    html.Div(id="save-figure-status", className="text-muted"),
                ]
            ),
        ],
        className="scb-viewlabel-card",
    )
