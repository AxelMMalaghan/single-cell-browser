from __future__ import annotations

import base64
import io
import json
import logging
import zipfile
from pathlib import Path
from typing import TYPE_CHECKING

import dash
import dash_bootstrap_components as dbc
from dash import ALL, Input, Output, State, dcc, html, no_update

from sc_browser.core.filter_state import FilterState
from sc_browser.metadata_io.model import session_from_dict, touch_session, session_to_dict
from sc_browser.services.session_service import delete_figure as svc_delete_figure
from sc_browser.ui.ids import IDs
from sc_browser.validation.errors import ValidationError
from sc_browser.validation.session_validation import validate_session_import_dict

if TYPE_CHECKING:
    from sc_browser.ui.config import AppConfig

logger = logging.getLogger(__name__)

MAX_IMPORTED_FIGURES = 100


def register_reports_callbacks(app: dash.Dash, ctx: AppConfig) -> None:
    """
    Callbacks for the Reports tab: render summary, list of saved figures,
    handle Delete, and import/export via shared session metadata.
    """

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
                            id={"type": IDs.Pattern.REPORTS_DELETE, "index": fig.id},
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

        session = svc_delete_figure(session, figure_id=figure_id)
        touch_session(session)
        return session_to_dict(session)

    # ---------------------------------------------------------
    # Export session: metadata JSON + all figure images as a ZIP
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

        session_id = session.session_id or active_session_id or "session"

        export_root = getattr(ctx.export_service, "output_root", None) or getattr(ctx.export_service, "_output_root", None)
        if export_root is None:
            logger.exception("ExportService has no output_root/_output_root")
            raise dash.exceptions.PreventUpdate

        export_root = Path(export_root)
        session_dir = export_root / session_id

        logger.info(
            "Export ZIP requested: session_id=%s export_root=%s session_dir=%s",
            session_id,
            export_root,
            session_dir,
        )

        def _write_zip(buf: io.BytesIO):
            with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
                meta_dict = session_to_dict(session)
                zf.writestr("metadata.json", json.dumps(meta_dict, indent=2).encode("utf-8"))

                if not session_dir.exists():
                    zf.writestr(
                        "FIGURES_MISSING.txt",
                        f"No figures directory found at: {session_dir}\n".encode("utf-8"),
                    )
                    logger.warning("No figures directory at %s (export will contain metadata only)", session_dir)
                    return

                pngs = list(session_dir.glob("*.png"))
                if not pngs:
                    zf.writestr(
                        "FIGURES_MISSING.txt",
                        f"No PNGs found in: {session_dir}\n".encode("utf-8"),
                    )
                    logger.warning("No PNGs in %s (export will contain metadata only)", session_dir)
                    return

                for png_path in pngs:
                    zf.write(png_path, arcname=f"figures/{png_path.name}")

        filename = f"{session_id}_report.zip"
        return dcc.send_bytes(_write_zip, filename)

    # ---------------------------------------------------------
    # Import helpers
    # ---------------------------------------------------------
    def _fs_json(fig) -> str:
        """
        Normalise a figure's filter_state to a stable JSON string so
        exact duplicates can be detected even if list ordering differs.
        """
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

    def _merge_and_normalise_figures(existing_figures, imported_figures):
        existing = list(existing_figures or [])
        imported = list(imported_figures or [])

        # 1) Drop only exact dupes (relative to existing + what we've already kept)
        seen = {
            (getattr(f, "dataset_key", None), getattr(f, "view_id", None), _fs_json(f), getattr(f, "label", None))
            for f in existing
        }
        merged = list(existing)

        for f in imported:
            key = (getattr(f, "dataset_key", None), getattr(f, "view_id", None), _fs_json(f), getattr(f, "label", None))
            if key in seen:
                continue
            seen.add(key)
            merged.append(f)

        # 2) PURELY sequential IDs, always (fig-0001..fig-N)
        for i, f in enumerate(merged, start=1):
            f.id = f"fig-{i:04d}"

        return merged

    # ---------------------------------------------------------
    # Import session metadata from metadata.json
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
        except ValueError:
            return (
                no_update,
                "Import failed: could not read the uploaded file. "
                "Please upload the metadata.json file exported from this app.",
            )

        try:
            raw_bytes = base64.b64decode(b64data)
            imported_dict = json.loads(raw_bytes.decode("utf-8"))
        except Exception as e:
            return (
                no_update,
                f"Import failed: invalid JSON ({e}). "
                "If you exported a ZIP, unzip it and upload metadata.json.",
            )

        try:
            validate_session_import_dict(imported_dict)
        except ValidationError as e:
            return (
                no_update,
                "Import failed:\n" + "\n".join(f"- {issue.message}" for issue in e.issues),
            )

        imported_session = session_from_dict(imported_dict)
        if imported_session is None:
            return (
                no_update,
                "Import failed: file didn’t look like session metadata. "
                "Upload the exported metadata.json.",
            )

        if len(imported_session.figures) > MAX_IMPORTED_FIGURES:
            logger.warning(
                "Truncating imported session figures: %d > MAX_IMPORTED_FIGURES=%d",
                len(imported_session.figures),
                MAX_IMPORTED_FIGURES,
            )
            imported_session.figures = imported_session.figures[:MAX_IMPORTED_FIGURES]

        existing_session = session_from_dict(session_data)
        existing_figures = existing_session.figures if existing_session else []

        merged_figures = _merge_and_normalise_figures(existing_figures, imported_session.figures)

        merged_session = existing_session or imported_session
        merged_session.figures = merged_figures
        touch_session(merged_session)

        # sanity log: should always be unique now
        ids = [f.id for f in merged_session.figures]
        logger.info("import: figures=%d unique=%d", len(ids), len(set(ids)))

        banner = (
            f"Imported {len(imported_session.figures)} figure(s). "
            f"Session now has {len(merged_figures)} figure(s)."
        )
        return session_to_dict(merged_session), banner
