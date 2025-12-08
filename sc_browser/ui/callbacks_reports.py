# sc_browser/ui/callbacks_reports.py

from __future__ import annotations

from typing import TYPE_CHECKING

import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, html, ALL

from sc_browser.export.model import session_from_dict, touch_session, session_to_dict

if TYPE_CHECKING:
    from .context import AppContext


def register_reports_callbacks(app: dash.Dash, ctx: "AppContext") -> None:
    """
    Callbacks for the Reports tab: render summary, list of saved figures,
    and handle Load/Delete actions.
    """

    # ---------------------------------------------------------
    # Render summary + table of figures
    # ---------------------------------------------------------
    @app.callback(
        Output("reports-summary-text", "children"),
        Output("reports-figure-list", "children"),
        Input("session-metadata", "data"),
    )
    def update_reports_list(session_data):
        session = session_from_dict(session_data)

        # No session or no figures yet
        if session is None or not session.figures:
            summary = "No figures saved in this session yet."
            empty = html.Div(
                "No figures saved yet. Go to the Explore tab and click 'Save figure'.",
                className="text-muted small mt-2",
            )
            return summary, empty

        n = len(session.figures)
        summary = f"{n} saved figure{'s' if n != 1 else ''} in this session."

        header = html.Thead(
            html.Tr(
                [
                    html.Th("ID"),
                    html.Th("Label"),
                    html.Th("Dataset"),
                    html.Th("View"),
                    html.Th("Created at"),
                    html.Th("Actions"),
                ]
            )
        )

        body_rows = []
        for fig in session.figures:
            body_rows.append(
                html.Tr(
                    [
                        html.Td(fig.id),
                        html.Td(fig.label or ""),
                        html.Td(fig.dataset_key),
                        html.Td(fig.view_id),
                        html.Td(fig.created_at),
                        html.Td(
                            [
                                dbc.Button(
                                    "Load",
                                    id={"type": "reports-load-figure", "index": fig.id},
                                    color="primary",
                                    size="sm",
                                    className="me-2",
                                ),
                                dbc.Button(
                                    "Delete",
                                    id={"type": "reports-delete-figure", "index": fig.id},
                                    color="danger",
                                    outline=True,
                                    size="sm",
                                ),
                            ]
                        ),
                    ]
                )
            )

        table = dbc.Table(
            [header, html.Tbody(body_rows)],
            bordered=True,
            hover=True,
            size="sm",
            className="mt-2",
        )

        return summary, table



    @app.callback(
        Output("session-metadata", "data"),
        Input({"type": "reports-delete-figure", "index": ALL}, "n_clicks"),
        State("session-metadata", "data"),
        prevent_initial_call=True,
    )
    def delete_figure(n_clicks_list, session_data):
        # If nothing was ever clicked, don't touch state
        if not n_clicks_list or all((n is None or n == 0) for n in n_clicks_list):
            raise dash.exceptions.PreventUpdate

        ctx_trigger = dash.callback_context
        if not ctx_trigger.triggered:
            raise dash.exceptions.PreventUpdate

        # Figure out which button fired
        prop_id = ctx_trigger.triggered[0]["prop_id"].split(".")[0]
        import json
        trigger_id = json.loads(prop_id)   # id={"type": "...", "index": "<figure_id>"}
        figure_id = trigger_id["index"]

        session = session_from_dict(session_data)
        if session is None:
            raise dash.exceptions.PreventUpdate

        # Filter out the deleted figure
        session.figures = [f for f in session.figures if f.id != figure_id]

        touch_session(session)
        return session_to_dict(session)




    @app.callback(
        Output("active-figure-id", "data"),
        Output("view-tabs", "value"),
        Output("cluster-select", "value"),
        Output("condition-select", "value"),
        Output("sample-select", "value"),
        Output("celltype-select", "value"),
        Output("gene-select", "value"),
        Output("embedding-select", "value"),
        Output("dim-x-select", "value"),
        Output("dim-y-select", "value"),
        Output("dim-z-select", "value"),
        Output("options-checklist", "value"),
        Input({"type": "reports-load-figure", "index": ALL}, "n_clicks"),
        State("session-metadata", "data"),
        State("dataset-select", "value"),
        prevent_initial_call=True,
    )
    def load_figure(
        n_clicks_list,
        session_data,
        current_dataset_name,
    ):
        # Nothing clicked -> do nothing
        if not n_clicks_list or all((n is None or n == 0) for n in n_clicks_list):
            raise dash.exceptions.PreventUpdate

        ctx_trigger = dash.callback_context
        if not ctx_trigger.triggered:
            raise dash.exceptions.PreventUpdate

        # Determine which Load button fired
        import json
        prop_id = ctx_trigger.triggered[0]["prop_id"].split(".")[0]
        trigger_id = json.loads(prop_id)   # {"type": "reports-load-figure", "index": "<figure_id>"}
        figure_id = trigger_id["index"]

        session = session_from_dict(session_data)
        if session is None:
            raise dash.exceptions.PreventUpdate

        # Find the matching figure
        meta = next((f for f in session.figures if f.id == figure_id), None)
        if meta is None:
            raise dash.exceptions.PreventUpdate

        # v1: only load if this figure belongs to the currently selected dataset
        # meta.dataset_key is what you saved (we used ds.name or ds.key)
        # current_dataset_name is the dataset-select value (name)
        # If you used a separate Dataset.key, resolve it via ctx.dataset_by_key:
        ds = ctx.dataset_by_key.get(meta.dataset_key)
        if ds is None or ds.name != current_dataset_name:
            # Option: you could set a banner here telling user to switch dataset
            raise dash.exceptions.PreventUpdate

        fs = meta.filter_state or {}

        clusters = fs.get("clusters", [])
        conditions = fs.get("conditions", [])
        samples = fs.get("samples", [])
        cell_types = fs.get("cell_types", [])
        genes = fs.get("genes", [])
        embedding = fs.get("embedding", None)
        dim_x = fs.get("dim_x", 0)
        dim_y = fs.get("dim_y", 1)
        dim_z = fs.get("dim_z", 2)
        split_by_condition = fs.get("split_by_condition", False)
        is_3d = fs.get("is_3d", False)

        options = []
        if split_by_condition:
            options.append("split_by_condition")
        if is_3d:
            options.append("is_3d")

        return (
            figure_id,      # active-figure-id
            meta.view_id,   # view-tabs
            clusters,
            conditions,
            samples,
            cell_types,
            genes,
            embedding,
            dim_x,
            dim_y,
            dim_z,
            options,
        )


