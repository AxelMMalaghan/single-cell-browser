from __future__ import annotations

import base64
import json
import logging
from typing import TYPE_CHECKING

import dash
from dash import ALL, Input, Output, State, dcc, no_update

from sc_browser.core.metadata_model import SessionMetadata, now_iso
from sc_browser.ui.ids import IDs
from sc_browser.ui.layout.build_reports_panel import (
    build_empty_figures_message,
    build_figures_table,
)

if TYPE_CHECKING:
    from sc_browser.ui.config import AppConfig

logger = logging.getLogger(__name__)

MAX_IMPORTED_FIGURES = 100


def register_reports_callbacks(app: dash.Dash, ctx: AppConfig) -> None:
    # ---------------------------------------------------------
    # 1) Render summary + table of figures
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.REPORTS_SUMMARY_TEXT, "children"),
        Output(IDs.Control.REPORTS_FIGURE_LIST, "children"),
        Input(IDs.Store.SESSION_META, "data"),
    )
    def update_reports_list(session_data):
        """
        Updates the UI list. session_data is now a dict
        representing the SessionMetadata container.
        """
        if not session_data or not isinstance(session_data, dict):
            return (
                "No saved figures in this session yet.",
                build_empty_figures_message(),
            )

        figures = session_data.get("figures", [])
        if not figures:
            summary = "No saved figures in this session yet."
            return summary, build_empty_figures_message()

        n = len(figures)
        summary = f"{n} saved figure{'s' if n != 1 else ''} in this session."

        # build_figures_table expects a list of figure dictionaries
        table = build_figures_table(figures)

        return summary, table

    # ---------------------------------------------------------
    # 2) Delete figure from session
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Store.SESSION_META, "data", allow_duplicate=True),
        Input({"type": IDs.Pattern.REPORTS_DELETE, "index": ALL}, "n_clicks"),
        State(IDs.Store.SESSION_META, "data"),
        prevent_initial_call=True,
    )
    def delete_figure(_n_clicks_list, session_data):
        logger.info(
            "delete_figure fired: triggered=%s value=%s",
            dash.ctx.triggered_id,
            dash.ctx.triggered[0].get("value") if dash.ctx.triggered else None,
        )

        # Only act on a real click (not re-render / remount)
        if not _n_clicks_list or max((n or 0) for n in _n_clicks_list) == 0:
            raise dash.exceptions.PreventUpdate

        triggered = dash.ctx.triggered_id
        if not triggered or not isinstance(triggered, dict):
            raise dash.exceptions.PreventUpdate

        figure_id = triggered.get("index")
        if not figure_id or not session_data:
            raise dash.exceptions.PreventUpdate

        session = SessionMetadata.from_dict(session_data)
        new_figures = [f for f in session.figures if f.id != figure_id]
        if len(new_figures) == len(session.figures):
            raise dash.exceptions.PreventUpdate

        session.figures = new_figures
        session.updated_at = now_iso()
        return session.to_dict()

    # ---------------------------------------------------------
    # 3) Export figures as ZIP (Metadata + Images)
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.REPORTS_DOWNLOAD_SESSION, "data"),
        Input(IDs.Control.REPORTS_EXPORT_BTN, "n_clicks"),
        State(IDs.Store.SESSION_META, "data"),
        prevent_initial_call=True,
    )
    def export_session_bundle(n_clicks, session_data):
        """
        Calls ExportService to render all figures and package them
        into a ZIP with the metadata.json.
        """
        if not n_clicks or not session_data:
            raise dash.exceptions.PreventUpdate

        try:
            # Trigger the stateless export service to generate ZIP bytes
            # Note: ExportService.create_session_zip should be refactored
            # to accept the session dictionary directly or a SessionMetadata object.
            zip_bytes = ctx.export_service.create_session_zip(session_data)

            filename = f"report_{session_data.get('session_id', 'export')}.zip"
            return dcc.send_bytes(zip_bytes, filename)
        except Exception as err:
            logger.exception("ZIP Export failed")
            # In a real app, you might want to return a notification here
            raise dash.exceptions.PreventUpdate from err

    # ---------------------------------------------------------
    # 4) Import figures from JSON
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Store.SESSION_META, "data", allow_duplicate=True),
        Output(IDs.Control.REPORTS_IMPORT_BANNER, "children"),
        Input(IDs.Control.REPORTS_UPLOAD, "contents"),
        prevent_initial_call=True,
    )
    def import_session_metadata(contents):
        if contents is None:
            raise dash.exceptions.PreventUpdate

        if isinstance(contents, list):
            contents = contents[0] if contents else None

        try:
            _header, b64data = contents.split(",", 1)
            raw_bytes = base64.b64decode(b64data)
            imported_json = json.loads(raw_bytes.decode("utf-8"))

            # Use SessionMetadata.from_dict to validate structure
            # and reconstitute FilterState objects
            new_session = SessionMetadata.from_dict(imported_json)

        except Exception as e:
            logger.error(f"Import failed: {e}")
            return no_update, f"Import failed: Invalid session metadata file ({e})."

        # Validate datasets exist in current app context
        unknown_datasets = set()
        for fig in new_session.figures:
            dataset_key = fig.dataset_key
            is_known = (dataset_key in ctx.dataset_by_key) or (
                dataset_key in ctx.dataset_by_name
            )
            if not is_known:
                unknown_datasets.add(dataset_key)

        # Truncate if too many figures
        if len(new_session.figures) > MAX_IMPORTED_FIGURES:
            logger.warning(
                f"Truncating import: {len(new_session.figures)} > {MAX_IMPORTED_FIGURES}"
            )
            new_session.figures = new_session.figures[:MAX_IMPORTED_FIGURES]

        banner_text = f"Imported session '{new_session.session_id}' with {len(new_session.figures)} figure(s)."

        if unknown_datasets:
            missing_str = ", ".join(sorted(unknown_datasets))
            banner_text += (
                f" Warning: some figures reference unknown datasets ({missing_str})."
            )

        return new_session.to_dict(), banner_text
