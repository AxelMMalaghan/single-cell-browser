from __future__ import annotations

import logging
from dataclasses import replace
from typing import TYPE_CHECKING

import dash
from dash import Input, Output, State, exceptions, dcc

from sc_browser.core.filter_state import FilterState
from sc_browser.metadata_io.metadata_model import (
    generate_session_id,
    session_from_dict,
    session_to_dict,
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
    # 1. Saved Figure Dropdown (UI State)
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.SAVED_FIGURE_SELECT, "options"),
        Output(IDs.Control.SAVED_FIGURE_SELECT, "value"),
        Input(IDs.Store.SESSION_META, "data"),
        State(IDs.Control.SAVED_FIGURE_SELECT, "value"),
    )
    def update_saved_figure_dropdown(session_data, current_value):
        """Updates options list when session data changes while preserving user selection."""
        options = [{"label": "New view (current filters)", "value": NEW_FIGURE_VALUE}]

        session = session_from_dict(session_data)
        if session is None or not session.figures:
            return options, NEW_FIGURE_VALUE

        ids: list[str] = []
        for fig in session.figures:
            display = fig.label or fig.view_id
            display = f"{display} ({fig.dataset_key})"
            options.append({"label": display, "value": fig.id})
            ids.append(fig.id)

        # Keep what the user selected if it's still valid; otherwise default to New.
        value = current_value if (current_value in ids or current_value == NEW_FIGURE_VALUE) else NEW_FIGURE_VALUE
        return options, value

    # ---------------------------------------------------------
    # 2. Sync Active Figure ID (Logic State)
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Store.ACTIVE_FIGURE_ID, "data", allow_duplicate=True),
        Input(IDs.Control.SAVED_FIGURE_SELECT, "value"),
        prevent_initial_call=True,
    )
    def sync_active_figure_id(selected):
        """
        Transitions the app between 'Edit Mode' (Overwrite) and 'New Mode' (Append)
        immediately when the dropdown selection changes.
        """
        return None if (not selected or selected == NEW_FIGURE_VALUE) else selected

    @app.callback(
        Output(IDs.Store.ACTIVE_FIGURE_ID, "data", allow_duplicate=True),
        Input(IDs.Control.DATASET_SELECT, "value"),
        prevent_initial_call=True,
    )
    def clear_active_figure_on_dataset_change(_ds_val):
        """Prevents cross-dataset overwrites by clearing edit mode on dataset switch."""
        return None

    # ---------------------------------------------------------
    # 3. Load Figure Logic
    # ---------------------------------------------------------
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
    def load_figure_from_dropdown(n_clicks, figure_id, session_data):
        if not n_clicks or not figure_id or figure_id == NEW_FIGURE_VALUE:
            return (dash.no_update,) * 14 + (None,)

        session = session_from_dict(session_data)
        meta = next((f for f in session.figures if f.id == figure_id), None) if session else None
        if not meta:
            raise exceptions.PreventUpdate

        ds = ctx.dataset_by_key.get(meta.dataset_key) or ctx.dataset_by_name.get(meta.dataset_key)
        if not ds:
            raise exceptions.PreventUpdate

        state = FilterState.from_dict(meta.filter_state or {})
        state = replace(state, dataset_name=ds.name, view_id=meta.view_id)

        options = []
        if state.split_by_condition: options.append(OPT_SPLIT_BY_CONDITION)
        if state.is_3d: options.append(OPT_IS_3D)

        return (
            state.dataset_name, state.view_id, state.clusters, state.conditions,
            state.samples, state.cell_types, state.genes, state.embedding,
            state.dim_x, state.dim_y, state.dim_z, options, state.color_scale,
            meta.label or "", figure_id
        )

    # ---------------------------------------------------------
    # 4. Save/Overwrite Figure Logic
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Store.SESSION_META, "data"),
        Output(IDs.Store.ACTIVE_SESSION_ID, "data"),
        Output(IDs.Control.SAVE_FIGURE_STATUS, "children"),
        Output(IDs.Store.ACTIVE_FIGURE_ID, "data", allow_duplicate=True),
        Input(IDs.Control.SAVE_FIGURE_BTN, "n_clicks"),
        State(IDs.Store.FILTER_STATE, "data"),
        State(IDs.Control.FIGURE_LABEL_INPUT, "value"),
        State(IDs.Store.SESSION_META, "data"),
        State(IDs.Store.ACTIVE_SESSION_ID, "data"),
        State(IDs.Store.ACTIVE_FIGURE_ID, "data"),
        State(IDs.Control.SAVED_FIGURE_SELECT, "value"),  # Use dropdown state as source of truth
        prevent_initial_call=True,
    )
    def save_current_figure(
            n_clicks, fs_data, figure_label, session_data, active_session_id, active_figure_id, current_selection
    ):
        if not n_clicks or not fs_data:
            raise exceptions.PreventUpdate

        # Robust Intent Check: If dropdown says "New", we FORCE a new figure regardless of store state.
        effective_id = None if current_selection == NEW_FIGURE_VALUE else active_figure_id

        try:
            state = FilterState.from_dict(fs_data)
            ds = ctx.dataset_by_name.get(state.dataset_name)
            if not ds: return dash.no_update, dash.no_update, "Dataset unavailable.", dash.no_update

            active_session_id = active_session_id or generate_session_id()
            session = ctx.session_service.ensure_session(
                session_from_dict(session_data),
                session_id=active_session_id
            )

            session, saved_id, is_overwrite = ctx.session_service.save_figure(
                session,
                active_figure_id=effective_id,
                dataset_key=getattr(ds, "key", ds.name),
                view_id=state.view_id,
                filter_state=state.to_dict(),
                label=figure_label,
            )

            status = f"{'Updated' if is_overwrite else 'Saved'} '{figure_label or state.view_id}' to report"
            # Keep ID if updating (stay in edit mode); clear if new (reset to create mode).
            next_active = saved_id if is_overwrite else None

            return session_to_dict(session), active_session_id, status, next_active

        except Exception as e:
            logger.error("Save failed: %s", str(e))
            return dash.no_update, dash.no_update, "Internal Save Error", dash.no_update

    # ---------------------------------------------------------
    # 5. Export Logic
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.DOWNLOAD_DATA, "data"),
        Output(IDs.Control.DOWNLOAD_DATA_BTN, "disabled"),
        Input(IDs.Control.DOWNLOAD_DATA_BTN, "n_clicks"),
        State(IDs.Store.FILTER_STATE, "data"),
        prevent_initial_call=True,
    )
    def download_current_data(n_clicks, fs_data):
        if not n_clicks or not fs_data:
            raise exceptions.PreventUpdate

        state = FilterState.from_dict(fs_data)
        ds = ctx.dataset_by_name.get(state.dataset_name)
        if not ds: raise exceptions.PreventUpdate

        view = ctx.registry.create(state.view_id, ds)
        data = view.compute_data(state)

        if hasattr(data, "to_csv"):
            return dcc.send_data_frame(data.to_csv, f"{state.view_id}_{ds.name}.csv", index=False), False

        return dash.no_update, False