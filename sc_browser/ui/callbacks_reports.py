from __future__ import annotations

import base64
import json
import logging
import dash
import dash_bootstrap_components as dbc

from dash import Input, Output, State, html, ALL, no_update
from sc_browser.metadata_io.model import session_from_dict, touch_session, session_to_dict
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .config import AppConfig

logger = logging.getLogger(__name__)

MAX_IMPORTED_FIGURES = 100


def register_reports_callbacks(app: dash.Dash, ctx: AppConfig) -> None:
    """
    Callbacks for the Reports tab: render summary, list of saved figures,
    and handle Load/Delete actions, plus integration with the Explore
    saved-figure dropdown via shared session metadata.
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

        body_rows = [
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
            for fig in session.figures
        ]

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
        trigger_id = json.loads(prop_id)  # {"type": "...", "index": "<figure_id>"}
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
            # No figures -> no download.
            raise dash.exceptions.PreventUpdate

        session_id = session.session_id or active_session_id or "session"

        export_root = getattr(ctx.export_service, "output_root", None)
        if export_root is None:
            export_root = getattr(ctx.export_service, "_output_root", None)

        if export_root is None:
            logger.error("ExportService has no output_root/_output_root; aborting export_session")
            raise dash.exceptions.PreventUpdate

        if not isinstance(export_root, Path):
            export_root = Path(export_root)

        session_dir = export_root / session_id

        def _write_zip(buf: io.BytesIO):
            with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
                # 1) metadata.json
                meta_dict = session_to_dict(session)
                meta_bytes = json.dumps(meta_dict, indent=2).encode("utf-8")
                zf.writestr("metadata.json", meta_bytes)

                # 2) PNG images
                if session_dir.exists():
                    for png_path in session_dir.glob("*.png"):
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
        The file should be the metadata.json we exported earlier.
        """
        if contents is None:
            raise dash.exceptions.PreventUpdate

        try:
            _header, b64data = contents.split(",", 1)
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

        imported_session = session_from_dict(imported_dict)
        if imported_session is None:
            return (
                no_update,
                "Import failed: JSON did not look like session metadata. "
                "Make sure you are uploading the metadata.json file created by the export button.",
            )

        # Hard cap on figures
        n_total = len(imported_session.figures)
        if n_total > MAX_IMPORTED_FIGURES:
            logger.warning(
                "Truncating imported session figures: %d > MAX_IMPORTED_FIGURES=%d",
                n_total,
                MAX_IMPORTED_FIGURES,
            )
            imported_session.figures = imported_session.figures[:MAX_IMPORTED_FIGURES]

        existing_session = session_from_dict(session_data)

        def _dedupe_figures(figures):
            seen_ids = set()
            unique = []
            for fig in figures:
                fig_id = getattr(fig, "id", None) or getattr(fig, "figure_id", None)
                if fig_id is None or fig_id in seen_ids:
                    continue
                seen_ids.add(fig_id)
                unique.append(fig)
            return unique

        if existing_session is None:
            imported_session.figures = _dedupe_figures(imported_session.figures)
            merged = imported_session
            merged_from_existing = False
        else:
            merged = existing_session
            merged.figures = _dedupe_figures(existing_session.figures + imported_session.figures)
            merged_from_existing = True

        touch_session(merged)

        n_imported_effective = len(imported_session.figures)

        if n_imported_effective == 0:
            banner = (
                "Imported metadata.json, but it did not contain any saved figures. "
                "Check that you selected the correct file."
            )
        else:
            truncated_note = ""
            if n_total > MAX_IMPORTED_FIGURES:
                truncated_note = f" (truncated to {MAX_IMPORTED_FIGURES} to keep the session manageable)"

            if merged_from_existing:
                banner = (
                    f"Imported {n_imported_effective} figure"
                    f"{'s' if n_imported_effective != 1 else ''}"
                    f"{truncated_note} and merged them into the current session. "
                    "They are now available in the Reports table and the Saved figure dropdown."
                )
            else:
                banner = (
                    f"Imported {n_imported_effective} figure"
                    f"{'s' if n_imported_effective != 1 else ''}"
                    f"{truncated_note} into a new session. "
                    "You can now load them from the Reports tab or the Saved figure dropdown."
                )

        return session_to_dict(merged), banner
