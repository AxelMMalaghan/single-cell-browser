from __future__ import annotations

import base64
import json
from typing import TYPE_CHECKING

import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, html, ALL, no_update

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
        Output("session-metadata", "data", allow_duplicate=True),
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
        Output("cluster-select", "value", allow_duplicate=True),
        Output("condition-select", "value", allow_duplicate=True),
        Output("sample-select", "value", allow_duplicate=True),
        Output("celltype-select", "value", allow_duplicate=True),
        Output("gene-select", "value", allow_duplicate=True),
        Output("embedding-select", "value", allow_duplicate=True),
        Output("dim-x-select", "value", allow_duplicate=True),
        Output("dim-y-select", "value", allow_duplicate=True),
        Output("dim-z-select", "value", allow_duplicate=True),
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


        ds = ctx.dataset_by_key.get(meta.dataset_key)
        if ds is None or ds.name != current_dataset_name:

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

    # ---------------------------------------------------------
    # Export session: metadata JSON + all figure images as a ZIP
    # ---------------------------------------------------------
    @app.callback(
        Output("reports-download-session", "data"),
        Input("reports-export-btn", "n_clicks"),
        State("session-metadata", "data"),
        State("active-session-id", "data"),
        prevent_initial_call=True,
    )
    def export_session(n_clicks, session_data, active_session_id):
        from dash import dcc as dash_dcc
        import io
        import json
        import zipfile
        from pathlib import Path

        if not n_clicks:
            raise dash.exceptions.PreventUpdate

        session = session_from_dict(session_data)
        if session is None or not session.figures:
            # nothing to export
            raise dash.exceptions.PreventUpdate

        # Use the real session_id from metadata; fall back to active_session_id if needed
        session_id = session.session_id or active_session_id or "session"

        # Root where images were written by GraphExportService
        # We set this as ctx.export_service._output_root in the service __init__
        export_root: Path = ctx.export_service._output_root  # yes, accessing "private" field
        session_dir = export_root / session_id

        # Build ZIP in memory
        def _write_zip(buf: io.BytesIO):
            with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
                # 1) metadata.json
                from sc_browser.export.model import session_to_dict  # adjust import if needed

                meta_dict = session_to_dict(session)
                meta_bytes = json.dumps(meta_dict, indent=2).encode("utf-8")
                zf.writestr("metadata.json", meta_bytes)

                # 2) all PNG images for this session, if the directory exists
                if session_dir.exists():
                    for png_path in session_dir.glob("*.png"):
                        # store under `figures/<filename>`
                        arcname = f"figures/{png_path.name}"
                        zf.write(png_path, arcname=arcname)

        filename = f"{session_id}_report.zip"
        return dash_dcc.send_bytes(_write_zip, filename)


    # ---------------------------------------------------------
    # Import session metadata from metadata.json
    # ---------------------------------------------------------
    @app.callback(
        Output("session-metadata", "data", allow_duplicate=True),
        Output("reports-import-banner", "children"),
        Input("reports-upload", "contents"),
        State("session-metadata", "data"),
        prevent_initial_call=True,
    )
    def import_session_metadata(contents, session_data):
        """
        Accepts a single JSON file from the Reports upload component.
        The file should be the metadata.json we exported earlier
        (i.e. output of session_to_dict).
        """
        if contents is None:
            raise dash.exceptions.PreventUpdate

        # contents is "data:application/json;base64,XXXXX"
        try:
            header, b64data = contents.split(",", 1)
        except ValueError:
            return no_update, "Import failed: could not parse upload payload."

        try:
            raw_bytes = base64.b64decode(b64data)
            raw_text = raw_bytes.decode("utf-8")
            imported_dict = json.loads(raw_text)
        except Exception as e:
            return no_update, f"Import failed: invalid JSON ({e})."

        # Rebuild SessionMetadata from uploaded JSON
        imported_session = session_from_dict(imported_dict)
        if imported_session is None:
            return no_update, "Import failed: JSON did not look like SessionMetadata."

        n_imported = len(imported_session.figures)

        # Merge into existing session if present; otherwise use imported as-is
        existing_session = session_from_dict(session_data)

        if existing_session is None:
            merged = imported_session
        else:
            # simple append; you can dedupe by id if you want
            existing_session.figures.extend(imported_session.figures)
            merged = existing_session

        touch_session(merged)
        banner = f"Imported {n_imported} figure{'s' if n_imported != 1 else ''} from metadata.json."

        return session_to_dict(merged), banner
