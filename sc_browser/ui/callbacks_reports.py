from __future__ import annotations

import base64
import json
from typing import TYPE_CHECKING

import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, html, ALL, no_update

from sc_browser.metadata_io.model import session_from_dict, touch_session, session_to_dict

if TYPE_CHECKING:
    from .context import AppContext


def register_reports_callbacks(app: dash.Dash, ctx: "AppContext") -> None:
    """
    Callbacks for the Reports tab: render summary, list of saved figures,
    and handle Load/Delete actions, plus integration with the Explore
    "Saved figure" dropdown.
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
            summary = "No saved figures in this session yet."
            empty = html.Div(
                [
                    html.Div(
                        "You haven’t saved any figures yet.",
                        className="fw-semibold",
                    ),
                    html.Div(
                        "Go to the Explore tab, configure a plot, and click “Save figure” to start building a report.",
                        className="text-muted small mt-1",
                    ),
                ],
                className="mt-2",
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
                        html.Td(fig.label or html.Span("No label", className="text-muted")),
                        html.Td(fig.dataset_key),
                        html.Td(fig.view_id),
                        html.Td(fig.created_at or ""),
                        html.Td(
                            dbc.Button(
                                "Delete",
                                id={"type": "reports-delete-figure", "index": fig.id},
                                color="danger",
                                outline=True,
                                size="sm",
                            )
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

    # ---------------------------------------------------------
    # Delete figure from session
    # ---------------------------------------------------------
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
        trigger_id = json.loads(prop_id)  # id={"type": "...", "index": "<figure_id>"}
        figure_id = trigger_id["index"]

        session = session_from_dict(session_data)
        if session is None:
            raise dash.exceptions.PreventUpdate

        # Filter out the deleted figure
        session.figures = [f for f in session.figures if f.id != figure_id]

        touch_session(session)
        return session_to_dict(session)

    # ---------------------------------------------------------
    # Export session: metadata JSON + all figure images as a ZIP
    # ---------------------------------------------------------
    @app.callback(
        Output("reports-download-session", "data"),
        Input("reports-metadata_io-btn", "n_clicks"),
        State("session-metadata", "data"),
        State("active-session-id", "data"),
        prevent_initial_call=True,
    )
    def export_session(n_clicks, session_data, active_session_id):
        from dash import dcc as dash_dcc
        import io
        import zipfile
        from pathlib import Path

        if not n_clicks:
            raise dash.exceptions.PreventUpdate

        session = session_from_dict(session_data)
        if session is None or not session.figures:
            # No figures -> no download. Button will appear to do nothing,
            # but we avoid sending an empty/broken file.
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
            return (
                no_update,
                "Import failed: could not read the uploaded file. "
                "Please select the metadata.json file that was exported from this app.",
            )

        try:
            raw_bytes = base64.b64decode(b64data)
            raw_text = raw_bytes.decode("utf-8")
            imported_dict = json.loads(raw_text)
        except Exception as e:
            return (
                no_update,
                f"Import failed: the file is not valid JSON ({e}). "
                "If you exported a ZIP, unzip it first and upload the metadata.json file inside.",
            )

        # Rebuild SessionMetadata from uploaded JSON
        imported_session = session_from_dict(imported_dict)
        if imported_session is None:
            return (
                no_update,
                "Import failed: JSON did not look like session metadata. "
                "Make sure you are uploading the metadata.json file created by the export button.",
            )

        n_imported = len(imported_session.figures)

        # Merge into existing session if present; otherwise use imported as-is
        existing_session = session_from_dict(session_data)

        if existing_session is None:
            merged = imported_session
            merged_from_existing = False
        else:
            # simple append; you can dedupe by id if you want
            existing_session.figures.extend(imported_session.figures)
            merged = existing_session
            merged_from_existing = True

        touch_session(merged)

        if n_imported == 0:
            banner = (
                "Imported metadata.json, but it did not contain any saved figures. "
                "Check that you selected the correct file."
            )
        else:
            if merged_from_existing:
                banner = (
                    f"Imported {n_imported} figure{'s' if n_imported != 1 else ''} "
                    "and merged them into the current session. "
                    "They are now available in the Reports table and the Saved figure dropdown."
                )
            else:
                banner = (
                    f"Imported {n_imported} figure{'s' if n_imported != 1 else ''} "
                    "into a new session. "
                    "You can now load them from the Reports tab or the Saved figure dropdown."
                )

        return session_to_dict(merged), banner

    # ---------------------------------------------------------
    # Populate "Select View" dropdown in the Explore view panel
    # ---------------------------------------------------------
    @app.callback(
        Output("saved-figure-select", "options"),
        Output("saved-figure-select", "value"),
        Input("session-metadata", "data"),
        State("saved-figure-select", "value"),
    )
    def update_saved_figure_dropdown(session_data, current_value):
        session = session_from_dict(session_data)

        # Base option: "New view"
        options = [
            {"label": "New view (start from current filters)", "value": "__new__"},
        ]

        if session is None or not session.figures:
            # No figures: only "New view"
            return options, "__new__"

        # Add each saved figure
        for fig in session.figures:
            base_label = fig.label or fig.id
            pretty = f"{base_label}  –  {fig.view_id} · {fig.dataset_key}"
            options.append({"label": pretty, "value": fig.id})

        valid_values = {opt["value"] for opt in options}

        # Keep current selection if still valid, else default to "New view"
        if current_value in valid_values:
            value = current_value
        else:
            value = "__new__"

        return options, value

    # ---------------------------------------------------------
    # Load a saved figure directly from the "Select View" dropdown
    # ---------------------------------------------------------
    @app.callback(
        Output("dataset-select", "value", allow_duplicate=True),
        Output("active-figure-id", "data", allow_duplicate=True),
        Output("view-select", "value", allow_duplicate=True),
        Output("figure-label-input", "value", allow_duplicate=True),
        Output("cluster-select", "value", allow_duplicate=True),
        Output("condition-select", "value", allow_duplicate=True),
        Output("sample-select", "value", allow_duplicate=True),
        Output("celltype-select", "value", allow_duplicate=True),
        Output("gene-select", "value", allow_duplicate=True),
        Output("embedding-select", "value", allow_duplicate=True),
        Output("dim-x-select", "value", allow_duplicate=True),
        Output("dim-y-select", "value", allow_duplicate=True),
        Output("dim-z-select", "value", allow_duplicate=True),
        Output("options-checklist", "value", allow_duplicate=True),
        Input("saved-figure-select", "value"),
        State("session-metadata", "data"),
        prevent_initial_call=True,
    )
    def load_figure_from_dropdown(
        figure_id,
        session_data,
    ):
        # "New view" sentinel or nothing selected -> don't change anything
        if not figure_id or figure_id == "__new__":
            raise dash.exceptions.PreventUpdate

        session = session_from_dict(session_data)
        if session is None:
            raise dash.exceptions.PreventUpdate

        # Find the matching saved figure
        meta = next((f for f in session.figures if f.id == figure_id), None)
        if meta is None:
            raise dash.exceptions.PreventUpdate

        # Get the dataset this figure belongs to
        ds = ctx.dataset_by_key.get(meta.dataset_key)
        if ds is None:
            # Dataset not available in current app state
            raise dash.exceptions.PreventUpdate

        fs = meta.filter_state or {}

        clusters = fs.get("clusters", [])
        conditions = fs.get("conditions", [])
        samples = fs.get("samples", [])
        cell_types = fs.get("cell_types", [])
        genes = fs.get("genes", [])

        # ---- NORMALISE EMBEDDING KEY ----
        embedding = fs.get("embedding", None)
        if not embedding or embedding not in ds.adata.obsm:
            default_emb = getattr(ds, "embedding_key", None)
            if default_emb and default_emb in ds.adata.obsm:
                embedding = default_emb
            else:
                embedding = None

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

        # Force dataset-select to match the figure's dataset
        dataset_name = ds.name

        return (
            dataset_name,       # dataset-select.value
            figure_id,          # active-figure-id
            meta.view_id,       # view-select (View type)
            meta.label or "",   # Figure label
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
