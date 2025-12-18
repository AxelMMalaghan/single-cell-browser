from __future__ import annotations

import logging
from dataclasses import replace
from typing import TYPE_CHECKING, Any

import dash
import pandas as pd
from dash import Input, Output, State, exceptions, dcc, no_update

from sc_browser.core.filter_state import FilterState
from sc_browser.metadata_io.metadata_model import (
    FigureMetadata,
    generate_figure_id,
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
    from dash import no_update

    @app.callback(
        Output(IDs.Control.SAVED_FIGURE_SELECT, "options"),
        Output(IDs.Control.SAVED_FIGURE_SELECT, "value"),
        Input(IDs.Store.SESSION_META, "data"),
        Input(IDs.Control.DATASET_SELECT, "value"),
        State(IDs.Store.ACTIVE_FIGURE_ID, "data"),
    )
    def update_saved_figure_dropdown(figures_data, dataset_value, active_figure_id):
        """
        Update dropdown options and value based on saved figures.

        figures_data is a simple list of figure dicts.
        Resets to "New view" when dataset changes.
        """
        figures = figures_data if isinstance(figures_data, list) else []
        triggered = dash.ctx.triggered_id
        logger.info(f"DROPDOWN: {len(figures)} figures, active_figure_id={active_figure_id}, triggered={triggered}")

        options = [{"label": "New view (current filters)", "value": NEW_FIGURE_VALUE}]

        # Reset to "New view" when dataset changes
        if triggered == IDs.Control.DATASET_SELECT:
            # Build options but reset selection
            for fig in figures:
                if isinstance(fig, dict):
                    fig_id = fig.get("id", "")
                    label = fig.get("label") or fig.get("view_id", "unknown")
                    dataset = fig.get("dataset_key", "")
                    options.append({"label": f"{label} ({dataset})", "value": fig_id})
            return options, NEW_FIGURE_VALUE

        if not figures:
            return options, NEW_FIGURE_VALUE

        ids: set[str] = set()
        for fig in figures:
            if isinstance(fig, dict):
                fig_id = fig.get("id", "")
                label = fig.get("label") or fig.get("view_id", "unknown")
                dataset = fig.get("dataset_key", "")
                display = f"{label} ({dataset})"
                options.append({"label": display, "value": fig_id})
                ids.add(fig_id)

        # Use ACTIVE_FIGURE_ID as the source of truth
        if active_figure_id and active_figure_id in ids:
            return options, active_figure_id

        return options, NEW_FIGURE_VALUE

    # ---------------------------------------------------------
    # 2) Clear ACTIVE_FIGURE_ID when user selects "New view" from dropdown
    # ---------------------------------------------------------
    # Note: We intentionally do NOT sync ACTIVE_FIGURE_ID when user selects
    # an existing figure from the dropdown. The user must click "Load" to
    # actually load the figure. This prevents circular callback issues.
    @app.callback(
        Output(IDs.Store.ACTIVE_FIGURE_ID, "data", allow_duplicate=True),
        Input(IDs.Control.SAVED_FIGURE_SELECT, "value"),
        State(IDs.Store.ACTIVE_FIGURE_ID, "data"),
        prevent_initial_call=True,
    )
    def clear_active_figure_on_new_selection(selected: str | None, current_active: str | None):
        # Only clear ACTIVE_FIGURE_ID if user explicitly selects "New view"
        # Don't change it when selecting an existing figure (user must click Load)
        if selected == NEW_FIGURE_VALUE and current_active is not None:
            return None
        # For all other cases, don't change ACTIVE_FIGURE_ID
        raise exceptions.PreventUpdate

    @app.callback(
        Output(IDs.Store.ACTIVE_FIGURE_ID, "data", allow_duplicate=True),
        Input(IDs.Control.DATASET_SELECT, "value"),
        prevent_initial_call=True,
    )
    def clear_active_figure_on_dataset_change(_ds_val: str | None):
        return None

    # ---------------------------------------------------------
    # 3) Load fig from dropdown
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
    def load_figure_from_dropdown(n_clicks: int | None, figure_id: str | None, figures_data: list | None):
        if not n_clicks or not figure_id or figure_id == NEW_FIGURE_VALUE:
            return (no_update,) * 14 + (None,)

        figures = figures_data if isinstance(figures_data, list) else []
        meta = next((f for f in figures if isinstance(f, dict) and f.get("id") == figure_id), None)
        if not meta:
            raise exceptions.PreventUpdate

        dataset_key = meta.get("dataset_key", "")
        ds = ctx.dataset_by_key.get(dataset_key) or ctx.dataset_by_name.get(dataset_key)
        if not ds:
            raise exceptions.PreventUpdate

        filter_state = meta.get("filter_state") or {}
        state = FilterState.from_dict(filter_state)
        state = replace(state, dataset_name=ds.name, view_id=meta.get("view_id", ""))

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
            meta.get("label") or "",
            figure_id,
        )

    # ---------------------------------------------------------
    # 4) Save fig - stores a simple list of figure dicts
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Store.SESSION_META, "data", allow_duplicate=True),
        Output(IDs.Control.SAVE_FIGURE_STATUS, "children"),
        Output(IDs.Store.ACTIVE_FIGURE_ID, "data", allow_duplicate=True),
        Input(IDs.Control.SAVE_FIGURE_BTN, "n_clicks"),
        State(IDs.Store.FILTER_STATE, "data"),
        State(IDs.Control.FIGURE_LABEL_INPUT, "value"),
        State(IDs.Store.SESSION_META, "data"),
        State(IDs.Store.ACTIVE_FIGURE_ID, "data"),
        prevent_initial_call=True,
    )
    def save_current_figure(
            n_clicks: int | None,
            fs_data: dict[str, Any] | None,
            figure_label: str | None,
            figures_data: list | None,
            active_figure_id: str | None,
    ):
        if not n_clicks or not fs_data:
            raise exceptions.PreventUpdate

        try:
            state = FilterState.from_dict(fs_data)
            ds = ctx.dataset_by_name.get(state.dataset_name)
            if not ds:
                return no_update, "Dataset unavailable.", no_update

            # Simple list of figures - no session wrapper
            figures = list(figures_data) if isinstance(figures_data, list) else []
            logger.info(f"SAVE INPUT: {len(figures)} figures")

            # Determine if we're overwriting an existing figure
            existing_idx = None
            if active_figure_id:
                existing_idx = next(
                    (i for i, f in enumerate(figures) if isinstance(f, dict) and f.get("id") == active_figure_id),
                    None
                )
            is_overwrite = existing_idx is not None

            # Create the figure as a simple dict
            label_clean = str(figure_label).strip() if figure_label else None
            new_id = active_figure_id if is_overwrite else generate_figure_id()

            fig_dict = {
                "id": new_id,
                "dataset_key": getattr(ds, "key", ds.name),
                "view_id": state.view_id,
                "filter_state": state.to_dict(),
                "label": label_clean,
            }

            # Update or append
            if is_overwrite:
                figures[existing_idx] = fig_dict
            else:
                figures.append(fig_dict)

            logger.info(f"SAVE OUTPUT: {len(figures)} figures, is_overwrite={is_overwrite}")

            status = f"{'Updated' if is_overwrite else 'Saved'} '{label_clean or state.view_id}' to report"

            if is_overwrite:
                return figures, status, new_id

            return figures, status, None

        except Exception as e:
            logger.exception("Save failed")
            return no_update, f"Save Error: {e}", no_update

    # ---------------------------------------------------------
    # 5) Export Logic
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
            return dcc.send_data_frame(data.to_csv, f"{state.view_id}_{ds.name}.csv", index=False)

        if isinstance(data, dict) and state.view_id == "dataset_summary":
            c_counts = data.get("cluster_counts", pd.DataFrame())
            cond_counts = data.get("condition_counts", pd.DataFrame())

            combined = pd.concat(
                [
                    c_counts.assign(type="cluster").rename(columns={"cluster": "label"}),
                    cond_counts.assign(type="condition").rename(columns={"condition": "label"}),
                ]
            )
            return dcc.send_data_frame(combined.to_csv, f"summary_{ds.name}.csv", index=False)

        logger.warning("Unsupported download type for view=%s type=%s", state.view_id, type(data))
        raise exceptions.PreventUpdate
