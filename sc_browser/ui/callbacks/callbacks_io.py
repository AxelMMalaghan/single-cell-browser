from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

import dash
import pandas as pd
from dash import Input, Output, State, exceptions, dcc, no_update

from sc_browser.core.filter_state import FilterState
from sc_browser.core.metadata_model import (
    FigureMetadata,
    SessionMetadata,
    new_session_metadata,
    generate_figure_id,
    now_iso,
)
from sc_browser.ui.ids import IDs

if TYPE_CHECKING:
    from sc_browser.ui.config import AppConfig

logger = logging.getLogger(__name__)

NEW_FIGURE_VALUE = "__new__"
OPT_SPLIT_BY_CONDITION = "split_by_condition"
OPT_IS_3D = "is_3d"


def register_io_callbacks(app: dash.Dash, ctx: AppConfig) -> None:
    # ---------------------------------------------------------
    # 1) Saved-figure dropdown options/value
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.SAVED_FIGURE_SELECT, "options"),
        Output(IDs.Control.SAVED_FIGURE_SELECT, "value"),
        Output(IDs.Store.SESSION_META, "data", allow_duplicate=True),
        Input(IDs.Store.SESSION_META, "data"),
        Input(IDs.Control.DATASET_SELECT, "value"),
        State(IDs.Store.ACTIVE_FIGURE_ID, "data"),
        prevent_initial_call="initial_duplicate",
    )
    def update_saved_figure_dropdown(session_data, dataset_value, active_figure_id):
        """
        Update dropdown based on the structured SessionMetadata container.
        """
        session_update = no_update
        figures: list[FigureMetadata] = []
        if isinstance(session_data, dict):
            raw_figures = session_data.get("figures", [])
            needs_update = any(
                (not isinstance(fig, dict)) or (not fig.get("id"))
                for fig in raw_figures
            )
            session = SessionMetadata.from_dict(session_data)
            seen_ids: set[str] = set()
            for fig in session.figures:
                if fig.id in seen_ids:
                    fig.id = generate_figure_id()
                    needs_update = True
                seen_ids.add(fig.id)
            if needs_update:
                session.updated_at = now_iso()
                session_update = session.to_dict()
            figures = session.figures

        triggered = dash.ctx.triggered_id
        options = [{"label": "New view (current filters)", "value": NEW_FIGURE_VALUE}]

        # Build options for all figures in the session
        figure_ids = set()
        for fig in figures:
            fig_id = fig.id
            label = fig.label or fig.view_id or "unknown"
            dataset = fig.dataset_key or ""
            options.append({"label": f"{label} ({dataset})", "value": fig_id})
            figure_ids.add(fig_id)

        # Reset to "New view" when dataset changes
        if triggered == IDs.Control.DATASET_SELECT:
            return options, NEW_FIGURE_VALUE, session_update

        # Use ACTIVE_FIGURE_ID as source of truth if it exists in the current session
        if active_figure_id and active_figure_id in figure_ids:
            return options, active_figure_id, session_update

        return options, NEW_FIGURE_VALUE, session_update

    # ---------------------------------------------------------
    # 2) Manage ACTIVE_FIGURE_ID lifecycle
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Store.ACTIVE_FIGURE_ID, "data", allow_duplicate=True),
        Input(IDs.Control.SAVED_FIGURE_SELECT, "value"),
        State(IDs.Store.ACTIVE_FIGURE_ID, "data"),
        prevent_initial_call=True,
    )
    def clear_active_figure_on_new_selection(
        selected: str | None, current_active: str | None
    ):
        if selected == NEW_FIGURE_VALUE and current_active is not None:
            return None
        raise exceptions.PreventUpdate

    @app.callback(
        Output(IDs.Store.ACTIVE_FIGURE_ID, "data", allow_duplicate=True),
        Input(IDs.Control.DATASET_SELECT, "value"),
        prevent_initial_call=True,
    )
    def clear_active_figure_on_dataset_change(_ds_val: str | None):
        return None

    # ---------------------------------------------------------
    # 3) Load figure filters from the session store
    # ---------------------------------------------------------
    # NOTE: This callback only updates UI controls, NOT FILTER_STATE directly.
    # FILTER_STATE is derived by sync_filter_state_from_ui (single-writer pattern)
    # to avoid race conditions between callbacks.
    @app.callback(
        Output(IDs.Control.DATASET_SELECT, "value", allow_duplicate=True),
        Output(IDs.Control.VIEW_SELECT, "value", allow_duplicate=True),
        Output(IDs.Control.CLUSTER_SELECT, "value", allow_duplicate=True),
        Output(IDs.Control.CONDITION_SELECT, "value", allow_duplicate=True),
        Output(IDs.Control.SAMPLE_SELECT, "value", allow_duplicate=True),
        Output(IDs.Control.CELLTYPE_SELECT, "value", allow_duplicate=True),
        Output(IDs.Control.GENE_SELECT, "value", allow_duplicate=True),
        Output(IDs.Control.EMBEDDING_SELECT, "value", allow_duplicate=True),
        Output(IDs.Control.DIM_X, "value", allow_duplicate=True),
        Output(IDs.Control.DIM_Y, "value", allow_duplicate=True),
        Output(IDs.Control.DIM_Z_SELECT, "value", allow_duplicate=True),
        Output(IDs.Control.OPTIONS_CHECKLIST, "value", allow_duplicate=True),
        Output(IDs.Control.COLOUR_SCALE_SELECT, "value", allow_duplicate=True),
        Output(IDs.Control.FIGURE_LABEL_INPUT, "value", allow_duplicate=True),
        Output(IDs.Store.ACTIVE_FIGURE_ID, "data", allow_duplicate=True),
        Input(IDs.Control.SAVED_FIGURE_LOAD_BTN, "n_clicks"),
        State(IDs.Control.SAVED_FIGURE_SELECT, "value"),
        State(IDs.Store.SESSION_META, "data"),
        prevent_initial_call=True,
    )
    def load_figure_from_session(
        n_clicks: int | None, figure_id: str | None, session_data: dict | None
    ):
        if (
            not n_clicks
            or not figure_id
            or figure_id == NEW_FIGURE_VALUE
            or not session_data
        ):
            return (no_update,) * 14 + (None,)

        figures = session_data.get("figures", [])
        fig_dict = next((f for f in figures if f.get("id") == figure_id), None)
        if not fig_dict:
            raise exceptions.PreventUpdate

        # Reconstruct typed FigureMetadata to ensure FilterState logic works
        meta = FigureMetadata.from_dict(fig_dict)
        state = meta.filter_state

        options: list[str] = []
        if state.split_by_condition:
            options.append(OPT_SPLIT_BY_CONDITION)
        if state.is_3d:
            options.append(OPT_IS_3D)

        # Return UI control values only. sync_filter_state_from_ui will derive
        # FILTER_STATE from these values, ensuring single source of truth.
        return (
            state.dataset_name,
            state.view_id,
            state.clusters,
            state.conditions,
            state.samples,
            state.cell_types,
            state.genes,
            state.embedding,
            state.dim_x,
            state.dim_y,
            state.dim_z,
            options,
            state.color_scale,
            meta.label or "",
            figure_id,
        )

    # ---------------------------------------------------------
    # 4) Save Current Figure to SessionMetadata
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Store.SESSION_META, "data", allow_duplicate=True),
        Output(IDs.Control.SAVE_FIGURE_STATUS, "children"),
        Output(IDs.Store.ACTIVE_FIGURE_ID, "data", allow_duplicate=True),
        Input(IDs.Control.SAVE_FIGURE_BTN, "n_clicks"),
        State(IDs.Store.FILTER_STATE, "data"),
        State(IDs.Control.FIGURE_LABEL_INPUT, "value"),
        State(IDs.Store.SESSION_META, "data"),
        prevent_initial_call=True,
    )
    def save_current_figure(
        n_clicks: int | None,
        fs_data: dict[str, Any] | None,
        figure_label: str | None,
        session_data: dict | None,
    ):
        if not n_clicks or not fs_data:
            raise exceptions.PreventUpdate

        try:
            app_ver = getattr(ctx.global_config, "app_version", "0.0.0-dev")

            # FIX: Create a deep copy of session_data to avoid mutation issues
            # and properly handle the figures list
            if session_data:
                # Get existing figures list
                existing_figures = session_data.get("figures", [])
                # Reconstruct session from dict
                session = SessionMetadata.from_dict(session_data)
            else:
                # Create new session if none exists
                existing_figures = []
                session = new_session_metadata(app_version=app_ver)

            # 1. Always generate a fresh ID
            new_id = generate_figure_id()

            # 2. Reconstruct the state and label
            state = FilterState.from_dict(fs_data)
            label_clean = str(figure_label).strip() if figure_label else None

            # 3. Create the new figure metadata
            new_fig = FigureMetadata(
                id=new_id,
                dataset_key=state.dataset_name,
                view_id=state.view_id,
                filter_state=state,
                label=label_clean,
                created_at=now_iso(),
            )

            # 4. FIX: Explicitly create a new list with all existing figures plus the new one
            # This ensures we don't have reference issues
            session.figures = [
                FigureMetadata.from_dict(f) if isinstance(f, dict) else f
                for f in existing_figures
            ] + [new_fig]
            session.updated_at = now_iso()

            status = f"Saved '{label_clean or state.view_id}' to report"

            # 5. Return dash.no_update for ACTIVE_FIGURE_ID.
            # This keeps the UI dropdown on "New view", allowing the user
            # to click Save again immediately to append another figure.
            return session.to_dict(), status, None

        except Exception as e:
            logger.exception("Save failed")
            return no_update, f"Save Error: {e}", no_update

    # ---------------------------------------------------------
    # 5) Data Export (CSV)
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.DOWNLOAD_DATA, "data"),
        Input(IDs.Control.DOWNLOAD_DATA_BTN, "n_clicks"),
        State(IDs.Store.FILTER_STATE, "data"),
        prevent_initial_call=True,
    )
    def download_current_data(n_clicks: int | None, fs_data: dict[str, Any] | None):
        if not n_clicks or not fs_data:
            raise exceptions.PreventUpdate

        state = FilterState.from_dict(fs_data)
        ds = ctx.dataset_by_name.get(state.dataset_name)
        if not ds:
            raise exceptions.PreventUpdate

        view = ctx.registry.create(state.view_id, ds)
        data = view.compute_data(state)

        if isinstance(data, pd.DataFrame):
            return dcc.send_data_frame(
                data.to_csv, f"{state.view_id}_{ds.name}.csv", index=False
            )

        raise exceptions.PreventUpdate
