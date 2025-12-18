from __future__ import annotations

import base64
import json
import logging
from typing import TYPE_CHECKING

import dash
from dash import ALL, Input, Output, State, dcc, no_update

from sc_browser.ui.ids import IDs
# Import layout builders
from sc_browser.ui.layout.build_reports_panel import build_empty_figures_message, build_figures_table

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
    def update_reports_list(figures_data):
        figures = figures_data if isinstance(figures_data, list) else []

        if not figures:
            summary = "No saved figures in this session yet."
            return summary, build_empty_figures_message()

        n = len(figures)
        summary = f"{n} saved figure{'s' if n != 1 else ''} in this session."

        # Use the helper function to build the styled table
        table = build_figures_table(figures)

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
    def delete_figure(_n_clicks_list, figures_data):
        triggered = dash.ctx.triggered_id
        if not triggered or not isinstance(triggered, dict):
            raise dash.exceptions.PreventUpdate

        figure_id = triggered.get("index")
        if not figure_id:
            raise dash.exceptions.PreventUpdate

        figures = figures_data if isinstance(figures_data, list) else []
        if not figures:
            raise dash.exceptions.PreventUpdate

        # Remove the figure from the list
        new_figures = [f for f in figures if isinstance(f, dict) and f.get("id") != figure_id]

        return new_figures

    # ---------------------------------------------------------
    # Export figures as JSON
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.REPORTS_DOWNLOAD_SESSION, "data"),
        Input(IDs.Control.REPORTS_EXPORT_BTN, "n_clicks"),
        State(IDs.Store.SESSION_META, "data"),
        prevent_initial_call=True,
    )
    def export_session(n_clicks, figures_data):
        if not n_clicks:
            raise dash.exceptions.PreventUpdate

        figures = figures_data if isinstance(figures_data, list) else []
        if not figures:
            raise dash.exceptions.PreventUpdate

        # Export as simple JSON
        json_str = json.dumps(figures, indent=2)
        return dcc.send_string(json_str, "figures_export.json")

    # ---------------------------------------------------------
    # Import figures from JSON
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Store.SESSION_META, "data", allow_duplicate=True),
        Output(IDs.Control.REPORTS_IMPORT_BANNER, "children"),
        Input(IDs.Control.REPORTS_UPLOAD, "contents"),
        State(IDs.Store.SESSION_META, "data"),
        prevent_initial_call=True,
    )
    def import_session_metadata(contents, figures_data):
        if contents is None:
            raise dash.exceptions.PreventUpdate

        if isinstance(contents, list):
            contents = contents[0] if contents else None
            if contents is None:
                raise dash.exceptions.PreventUpdate

        try:
            _header, b64data = contents.split(",", 1)
            raw_bytes = base64.b64decode(b64data)
            imported = json.loads(raw_bytes.decode("utf-8"))
        except Exception as e:
            return (
                no_update,
                f"Import failed: could not read uploaded file ({e}).",
            )

        # Handle different import formats
        imported_figures = []
        if isinstance(imported, list):
            # Direct list of figures
            imported_figures = [f for f in imported if isinstance(f, dict)]
        elif isinstance(imported, dict):
            # Could be a session dict with "figures" key or a single figure
            if "figures" in imported:
                imported_figures = [f for f in imported.get("figures", []) if isinstance(f, dict)]
            elif "id" in imported and "dataset_key" in imported:
                # Single figure dict
                imported_figures = [imported]

        if not imported_figures:
            return (
                no_update,
                "Import failed: no valid figures found in file.",
            )

        # Validate datasets exist
        unknown_datasets = set()
        for fig in imported_figures:
            dataset_key = fig.get("dataset_key", "")
            is_known = (dataset_key in ctx.dataset_by_key) or (dataset_key in ctx.dataset_by_name)
            if not is_known:
                unknown_datasets.add(dataset_key)

        # Limit imports
        if len(imported_figures) > MAX_IMPORTED_FIGURES:
            logger.warning(f"Truncating imported figures: {len(imported_figures)} > {MAX_IMPORTED_FIGURES}")
            imported_figures = imported_figures[:MAX_IMPORTED_FIGURES]

        # Merge with existing figures
        existing_figures = list(figures_data) if isinstance(figures_data, list) else []
        merged_figures = existing_figures + imported_figures

        logger.info(f"Import successful: added {len(imported_figures)}, total {len(merged_figures)}")

        banner_text = f"Imported {len(imported_figures)} figure(s). Session now has {len(merged_figures)} figure(s)."

        if unknown_datasets:
            missing_str = ", ".join(sorted(unknown_datasets))
            banner_text += f" Warning: some figures reference unknown datasets ({missing_str})."

        return merged_figures, banner_text


