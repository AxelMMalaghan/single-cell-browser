from __future__ import annotations

import dash_bootstrap_components as dbc
from dash import dcc, html
from typing import List, Any

# We need IDs for the delete button in the table
from sc_browser.ui.ids import IDs


def build_reports_panel() -> dbc.Container:
    """
    Reports page:
    - Shows how many figures are in the current session
    - Lists all figures (populated via callbacks in reports-figure-list)
    - Allows metadata_io of session metadata (+ images)
    - Allows import of a metadata JSON
    """

    # -- 1. Action Buttons (Export / Import) --
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
                                    # This div is populated by callback
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

                    # Sub-header
                    html.Label("Saved figures", className="form-label"),

                    # The list of figures (table) goes here
                    html.Div(
                        id="reports-figure-list",
                        className="scb-reports-figure-list mt-2",
                        style={"overflowX": "auto"},  # Ensure scroll on small screens
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


def build_figures_table(figures: List[Any]) -> dbc.Table:
    """
    Builds a styled dbc.Table for the reports list.
    Mimics the look-and-feel of a Dash DataTable (system fonts, specific padding/colors)
    but keeps the flexibility of embedding dbc.Button components.
    """

    # Styles derived from your dash_table.DataTable snippet
    font_family = 'system-ui, -apple-system, BlinkMacSystemFont, "SF Pro Text", "Segoe UI", sans-serif'

    header_style = {
        "fontFamily": font_family,
        "fontSize": "12px",
        "fontWeight": "600",
        "backgroundColor": "#f3f4f6",
        "borderBottom": "1px solid #e5e7eb",
        "color": "#111827",
        "padding": "8px 12px",
        "whiteSpace": "nowrap",
    }

    cell_style = {
        "fontFamily": font_family,
        "fontSize": "12px",
        "padding": "6px 12px",
        "borderBottom": "1px solid #e5e7eb",
        "color": "#374151",
        "verticalAlign": "middle",
        "whiteSpace": "nowrap",
    }

    # -- Header --
    thead = html.Thead(
        html.Tr(
            [
                html.Th("ID", style=header_style),
                html.Th("Label", style=header_style),
                html.Th("Dataset", style=header_style),
                html.Th("View", style=header_style),
                html.Th("Created at", style=header_style),
                html.Th("Actions", style=header_style),
            ]
        )
    )

    # -- Body --
    rows = []
    for fig in figures:
        rows.append(
            html.Tr(
                [
                    html.Td(fig.id, style=cell_style),
                    html.Td(
                        fig.label or html.Span("No label", className="text-muted italic"),
                        style={**cell_style, "whiteSpace": "normal", "maxWidth": "200px"}
                    ),
                    html.Td(fig.dataset_key, style=cell_style),
                    html.Td(fig.view_id, style=cell_style),
                    html.Td(fig.created_at or "", style=cell_style),
                    html.Td(
                        dbc.Button(
                            "Delete",
                            id={"type": IDs.Pattern.REPORTS_DELETE, "index": fig.id},
                            color="danger",
                            outline=True,
                            size="sm",
                            style={
                                "fontSize": "11px",
                                "padding": "2px 8px",
                                "lineHeight": "1.2"
                            }
                        ),
                        style=cell_style,
                    ),
                ]
            )
        )

    return dbc.Table(
        [thead, html.Tbody(rows)],
        bordered=False,  # We use our own border styles
        hover=True,
        responsive=True,
        className="mb-0",  # Remove bottom margin
        style={"border": "1px solid #e5e7eb", "borderRadius": "4px"}
    )