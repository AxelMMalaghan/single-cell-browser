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
    # Saved figure dropdown (Authoritative)
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.SAVED_FIGURE_SELECT, "options"),
        Output(IDs.Control.SAVED_FIGURE_SELECT, "value"),
        Input(IDs.Store.SESSION_META, "data"),
        State(IDs.Control.SAVED_FIGURE_SELECT, "value"),
        State(IDs.Store.ACTIVE_FIGURE_ID, "data"),
    )
    def update_saved_figure_dropdown(session_data, current_value, active_figure_id):
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

        # Logic: Prioritize active figure ID if present in list
        if active_figure_id in ids:
            value = active_figure_id
        elif current_value in ids:
            value = current_value
        else:
            value = NEW_FIGURE_VALUE

        return options, value

    # ---------------------------------------------------------
    # Load Figure
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
        if not n_clicks:
            raise exceptions.PreventUpdate

        if not figure_id or figure_id == NEW_FIGURE_VALUE:
            return (dash.no_update,) * 14 + (None,)

        session = session_from_dict(session_data)
        if session is None:
            raise exceptions.PreventUpdate

        meta = next((f for f in session.figures if f.id == figure_id), None)
        if meta is None:
            raise exceptions.PreventUpdate

        ds = ctx.dataset_by_key.get(meta.dataset_key)
        # Fallback lookups
        if ds is None:
            ds = ctx.dataset_by_name.get(meta.dataset_key)
        if ds is None:
            ds = next((d for d in ctx.datasets if d.name == meta.dataset_key), None)

        if ds is None:
            logger.warning("Dataset not found for figure %s", figure_id)
            raise exceptions.PreventUpdate

        # Parse saved state
        try:
            state = FilterState.from_dict(meta.filter_state or {})
        except Exception:
            raise exceptions.PreventUpdate

        # Ensure embedding params are valid
        state = replace(state, dataset_name=ds.name, view_id=meta.view_id)
        if state.embedding and state.embedding not in ds.adata.obsm:
            default_emb = getattr(ds, "embedding_key", None)
            state = replace(state, embedding=default_emb if default_emb in ds.adata.obsm else None)

        # Recalc dims if embedding changed or invalid
        if state.embedding:
            try:
                labels = ds.get_embedding_labels(state.embedding)
            except KeyError:
                labels = []
            if labels:
                dx = state.dim_x if (state.dim_x is not None and state.dim_x < len(labels)) else 0
                dy = state.dim_y if (state.dim_y is not None and state.dim_y < len(labels)) else 1
                state = replace(state, dim_x=dx, dim_y=dy)

        # Build UI values
        options: list[str] = []
        if state.split_by_condition:
            options.append(OPT_SPLIT_BY_CONDITION)
        if state.is_3d:
            options.append(OPT_IS_3D)

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
    # Save current figure
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Store.SESSION_META, "data"),
        Output(IDs.Store.ACTIVE_SESSION_ID, "data"),
        Output(IDs.Control.SAVE_FIGURE_STATUS, "children"),
        Output(IDs.Store.ACTIVE_FIGURE_ID, "data"),
        Input(IDs.Control.SAVE_FIGURE_BTN, "n_clicks"),
        State(IDs.Store.FILTER_STATE, "data"),
        State(IDs.Control.FIGURE_LABEL_INPUT, "value"),
        State(IDs.Store.SESSION_META, "data"),
        State(IDs.Store.ACTIVE_SESSION_ID, "data"),
        State(IDs.Store.ACTIVE_FIGURE_ID, "data"),
        prevent_initial_call=True,
    )
    def save_current_figure(
            n_clicks, fs_data, figure_label, session_data, active_session_id, active_figure_id
    ):
        if not n_clicks:
            raise exceptions.PreventUpdate

        if not fs_data:
            return dash.no_update, dash.no_update, "No state to save.", dash.no_update

        try:
            state = FilterState.from_dict(fs_data)
        except Exception:
            return dash.no_update, dash.no_update, "Invalid state.", dash.no_update

        ds = ctx.dataset_by_name.get(state.dataset_name)
        if ds is None:
            return dash.no_update, dash.no_update, "Dataset unavailable.", dash.no_update

        if active_session_id is None:
            active_session_id = generate_session_id()

        # REFACTOR: Delegate session management to Service (loads from disk or creates new)
        session = ctx.session_service.ensure_session(
            session_from_dict(session_data),
            session_id=active_session_id,
        )

        # Use the key (stable ID) if available, else name
        ds_key = getattr(ds, "key", ds.name)

        # REFACTOR: Save via Service (automatically persists to disk)
        session, new_fig_id, is_overwrite = ctx.session_service.save_figure(
            session,
            active_figure_id=active_figure_id,
            dataset_key=ds_key,
            view_id=state.view_id,
            filter_state=state.to_dict(),
            label=figure_label,
        )

        status = f"Saved '{figure_label or state.view_id}'" if not is_overwrite else f"Updated '{figure_label or state.view_id}'"

        return session_to_dict(session), active_session_id, status, new_fig_id

    # ---------------------------------------------------------
    # Download CSV
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

        try:
            state = FilterState.from_dict(fs_data)
        except Exception:
            raise exceptions.PreventUpdate

        ds = ctx.dataset_by_name.get(state.dataset_name)
        if not ds:
            raise exceptions.PreventUpdate

        view = ctx.registry.create(state.view_id, ds)
        data = view.compute_data(state)

        # re-enable button (return False for disabled)
        if hasattr(data, "to_csv"):
            filename = f"{state.view_id}_{state.dataset_name}.csv"
            return dcc.send_data_frame(data.to_csv, filename, index=False), False

        return dash.no_update, False