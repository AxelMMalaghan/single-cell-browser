from __future__ import annotations

import base64
import json
import logging
from typing import TYPE_CHECKING

import dash
import dash_bootstrap_components as dbc
from dash import ALL, Input, Output, State, dcc, html, no_update

from sc_browser.core.filter_state import FilterState
from sc_browser.metadata_io.io import normalise_session_dict
from sc_browser.metadata_io.metadata_model import session_from_dict, touch_session, session_to_dict
from sc_browser.ui.ids import IDs
# Import the new table builder
from sc_browser.ui.layout.build_reports_panel import build_figures_table

if TYPE_CHECKING:
    from sc_browser.ui.config import AppConfig

logger = logging.getLogger(__name__)

MAX_IMPORTED_FIGURES = 100


def register_reports_callbacks(app: dash.Dash, ctx: AppConfig) -> None:
    # ---------------------------------------------------------
    # Render summary + table of figures
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.REPORTS_SUMMARY_TEXT, "children"),
        Output(IDs.Control.REPORTS_FIGURE_LIST, "children"),
        Input(IDs.Store.SESSION_META, "data"),
    )
    def update_reports_list(session_data):
        session = session_from_dict(session_data)

        if session is None or not session.figures:
            summary = "No saved figures in this session yet."
            empty = html.Div(
                [
                    html.Div("You haven’t saved any figures yet.", className="fw-semibold"),
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

        # Use the new helper function to build the styled table
        table = build_figures_table(session.figures)

        return summary, table

    # ---------------------------------------------------------
    # Delete figure from session
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Store.SESSION_META, "data", allow_duplicate=True),
        Input({"type": IDs.Pattern.REPORTS_DELETE, "index": ALL}, "n_clicks"),
        State(IDs.Store.SESSION_META, "data"),
        prevent_initial_call=True,
    )
    def delete_figure(_n_clicks_list, session_data):
        triggered = dash.ctx.triggered_id
        if not triggered or not isinstance(triggered, dict):
            raise dash.exceptions.PreventUpdate

        figure_id = triggered.get("index")
        if not figure_id:
            raise dash.exceptions.PreventUpdate

        session = session_from_dict(session_data)
        if session is None:
            raise dash.exceptions.PreventUpdate

        # REFACTOR: Use SessionService (persists deletion to disk)
        session = ctx.session_service.delete_figure(session, figure_id=figure_id)

        return session_to_dict(session)

    # ---------------------------------------------------------
    # Export session
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.REPORTS_DOWNLOAD_SESSION, "data"),
        Input(IDs.Control.REPORTS_EXPORT_BTN, "n_clicks"),
        State(IDs.Store.SESSION_META, "data"),
        State(IDs.Store.ACTIVE_SESSION_ID, "data"),
        prevent_initial_call=True,
    )
    def export_session(n_clicks, session_data, active_session_id):
        if not n_clicks:
            raise dash.exceptions.PreventUpdate

        session = session_from_dict(session_data)
        if session is None or not session.figures:
            raise dash.exceptions.PreventUpdate

        if not session.session_id and active_session_id:
            session.session_id = active_session_id

        # REFACTOR: Use ExportService (reads files from storage)
        try:
            zip_bytes = ctx.export_service.create_session_zip(session)
        except Exception:
            logger.exception("Failed to create session ZIP")
            raise dash.exceptions.PreventUpdate

        filename = f"{session.session_id}_report.zip"
        return dcc.send_bytes(zip_bytes, filename)

    # ---------------------------------------------------------
    # Import helpers
    # ---------------------------------------------------------
    def _fs_json(fig) -> str:
        fs = getattr(fig, "filter_state", None)
        if isinstance(fs, FilterState):
            fs = fs.to_dict()
        elif fs is None:
            fs = {}
        elif not isinstance(fs, dict):
            to_dict = getattr(fs, "to_dict", None)
            fs = to_dict() if callable(to_dict) else {}

        for k in ("clusters", "conditions", "samples", "cell_types", "genes"):
            v = fs.get(k)
            if isinstance(v, list):
                fs[k] = sorted(v)

        return json.dumps(fs, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


    # ---------------------------------------------------------
    # Import session metadata
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Store.SESSION_META, "data", allow_duplicate=True),
        Output(IDs.Control.REPORTS_IMPORT_BANNER, "children"),
        Input(IDs.Control.REPORTS_UPLOAD, "contents"),
        State(IDs.Store.SESSION_META, "data"),
        prevent_initial_call=True,
    )
    def import_session_metadata(contents, session_data):
        if contents is None:
            raise dash.exceptions.PreventUpdate

        if isinstance(contents, list):
            contents = contents[0] if contents else None
            if contents is None:
                raise dash.exceptions.PreventUpdate

        try:
            _header, b64data = contents.split(",", 1)
            raw_bytes = base64.b64decode(b64data)
            imported_dict = json.loads(raw_bytes.decode("utf-8"))
        except Exception as e:
            return (
                no_update,
                f"Import failed: could not read uploaded file ({e}).",
            )

        try:
            imported_session = normalise_session_dict(imported_dict)
        except Exception as e:
            logger.warning("Normalization failed for imported JSON: %s", e)
            return (
                no_update,
                "Import failed: file format not recognized (must be a valid session JSON or list of figures).",
            )

        # -------------------------------------------------------
        # VALIDATION: Check if datasets exist in current config
        # -------------------------------------------------------
        unknown_datasets = set()
        for fig in imported_session.figures:
            # Check key manager first, then name manager (fallback)
            # Checks existence in config without triggering a load
            is_known = (fig.dataset_key in ctx.dataset_by_key) or \
                       (fig.dataset_key in ctx.dataset_by_name)

            if not is_known:
                unknown_datasets.add(fig.dataset_key)

        if len(imported_session.figures) > MAX_IMPORTED_FIGURES:
            logger.warning(
                "Truncating imported figures: %d > MAX_IMPORTED_FIGURES=%d",
                len(imported_session.figures),
                MAX_IMPORTED_FIGURES,
            )
            imported_session.figures = imported_session.figures[:MAX_IMPORTED_FIGURES]

        existing_session = session_from_dict(session_data)
        existing_figures = existing_session.figures if existing_session else []



        merged_session = existing_session or imported_session

        touch_session(merged_session)

        # Persist change
        if ctx.session_service:
            ctx.session_service.persist_session(merged_session)

        ids = [f.id for f in merged_session.figures]
        logger.info("Import successful: merged total figures=%d", len(ids))

        banner_text = (
            f"Imported {len(imported_session.figures)} figure(s). "
            f"Session now has {len(existing_figures)} figure(s)."
        )

        # Append warning if orphans detected
        if unknown_datasets:
            missing_str = ", ".join(sorted(unknown_datasets))
            banner_text += f" Warning: {len(unknown_datasets)} figure(s) reference unknown datasets ({missing_str}) and may not load."

        # Use a simpler styling for the banner text if it contains a warning
        return session_to_dict(merged_session), banner_text


