from __future__ import annotations

import logging
from dataclasses import replace
from typing import Any, Optional, TYPE_CHECKING

import dash
import pandas as pd
import plotly.graph_objs as go
from dash import Input, Output, State, exceptions, dcc as dash_dcc, html, no_update

from sc_browser.core.filter_state import FilterState
from sc_browser.metadata_io.metadata_model import (
    generate_session_id,
    session_from_dict,
    session_to_dict,
)
from sc_browser.services.session_service import ensure_session, save_figure
from sc_browser.ui.helpers import get_filter_dropdown_options
from sc_browser.ui.ids import IDs

if TYPE_CHECKING:
    from sc_browser.ui.config import AppConfig

logger = logging.getLogger(__name__)

NEW_FIGURE_VALUE = "__new__"
OPT_SPLIT_BY_CONDITION = "split_by_condition"
OPT_IS_3D = "is_3d"


# -----------------------------------------------------------------------------
# Small utilities
# -----------------------------------------------------------------------------
def _message_figure(title: str, details: Optional[str] = None) -> go.Figure:
    fig = go.Figure()
    text = title if details is None else f"{title}\n\n{details}"
    fig.add_annotation(
        text=text,
        showarrow=False,
        xref="paper",
        yref="paper",
        x=0.5,
        y=0.5,
    )
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    fig.update_layout(margin=dict(l=40, r=40, t=40, b=40))
    return fig


def _error_figure(details: str) -> go.Figure:
    return _message_figure("Something went wrong while rendering this view.", details)


def _parse_filter_state(raw: Any) -> FilterState:
    """
    Strict parse: returns a FilterState or raises.
    Centralises the boundary conversion without "safe optional" behaviour.
    """
    if raw is None:
        raise ValueError("filter_state is None")
    if not isinstance(raw, dict):
        raise TypeError(f"filter_state must be dict, got {type(raw)}")
    return FilterState.from_dict(raw)


def _view_label(ctx: "AppConfig", view_id: str) -> str:
    if not view_id:
        return "No view"
    try:
        for cls in ctx.registry.all_classes():
            if getattr(cls, "id", None) == view_id:
                return getattr(cls, "label", view_id)
    except Exception:
        pass
    return view_id


def register_explore_callbacks(app: dash.Dash, ctx: "AppConfig") -> None:
    # ---------------------------------------------------------
    # Sidebar metadata
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.SIDEBAR_DATASET_NAME, "children"),
        Output(IDs.Control.SIDEBAR_DATASET_META, "children"),
        Input(IDs.Control.DATASET_SELECT, "value"),
    )
    def update_sidebar_dataset_summary(dataset_name: str | None):
        ds = ctx.dataset_by_name.get(dataset_name) if dataset_name else None
        if ds is None:
            return "No dataset", "0 cells · 0 genes"
        return ds.name, f"{ds.adata.n_obs} cells · {ds.adata.n_vars} genes"

    # ---------------------------------------------------------
    # Update dropdown options for filters (per dataset)
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.CLUSTER_SELECT, "options"),
        Output(IDs.Control.CONDITION_SELECT, "options"),
        Output(IDs.Control.SAMPLE_SELECT, "options"),
        Output(IDs.Control.CELLTYPE_SELECT, "options"),
        Output(IDs.Control.EMBEDDING_SELECT, "options"),
        Input(IDs.Control.DATASET_SELECT, "value"),
    )
    def update_filters(dataset_name: str | None):
        ds = ctx.dataset_by_name.get(dataset_name) if dataset_name else None
        if ds is None:
            empty: list[dict] = []
            return empty, empty, empty, empty, empty
        return get_filter_dropdown_options(ds)

    # ---------------------------------------------------------
    # Live gene search (prefix, keep selected) — uses Dataset.valid_sets()
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.GENE_SELECT, "options"),
        Input(IDs.Control.GENE_SELECT, "search_value"),
        Input(IDs.Control.DATASET_SELECT, "value"),
        State(IDs.Control.GENE_SELECT, "value"),
    )
    def update_gene_options(search_value, dataset_name, selected_genes):
        ds = ctx.dataset_by_name.get(dataset_name) if dataset_name else None
        if ds is None:
            return []

        valid = ds.valid_sets()
        all_genes = sorted(valid.genes)
        gene_set = valid.genes

        selected_genes = [g for g in (selected_genes or []) if g in gene_set]

        if not search_value or not str(search_value).strip():
            return [{"label": g, "value": g} for g in selected_genes]

        query = str(search_value).upper()
        matches = [g for g in all_genes if g.upper().startswith(query)]
        suggestions = matches[:50]

        union = list(dict.fromkeys(selected_genes + suggestions))
        return [{"label": g, "value": g} for g in union]

    # ---------------------------------------------------------
    # Hide/show filters based on active view + dataset capabilities
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.CLUSTER_FILTER_CONTAINER, "style"),
        Output(IDs.Control.CONDITION_FILTER_CONTAINER, "style"),
        Output(IDs.Control.SAMPLE_FILTER_CONTAINER, "style"),
        Output(IDs.Control.CELLTYPE_FILTER_CONTAINER, "style"),
        Output(IDs.Control.GENE_FILTER_CONTAINER, "style"),
        Output(IDs.Control.EMBEDDING_FILTER_CONTAINER, "style"),
        Output(IDs.Control.DIM_FILTER_CONTAINER, "style"),
        Output(IDs.Control.OPTIONS_CONTAINER, "style"),
        Output(IDs.Control.OPTIONS_CHECKLIST, "options"),
        Output(IDs.Control.COLOUR_SCALE_SELECT, "style"),
        Input(IDs.Control.VIEW_SELECT, "value"),
        State(IDs.Control.DATASET_SELECT, "value"),
    )
    def update_filter_visibility(view_id: str | None, dataset_name: str | None):
        def style(flag: bool) -> dict:
            return {} if flag else {"display": "none"}

        def default_options():
            return [
                {"label": "Split by condition", "value": OPT_SPLIT_BY_CONDITION},
                {"label": "3D view", "value": OPT_IS_3D},
            ]

        ds = ctx.dataset_by_name.get(dataset_name) if dataset_name else None
        if ds is None or not view_id:
            options = default_options()
            return (
                style(True),
                style(True),
                style(True),
                style(True),
                style(True),
                style(True),
                style(True),
                style(True),
                options,
                style(False),
            )

        view = ctx.registry.create(view_id, ds)
        profile = getattr(view, "filter_profile", None)
        if profile is None:
            options = default_options()
            return (
                style(True),
                style(True),
                style(True),
                style(True),
                style(True),
                style(True),
                style(True),
                style(True),
                options,
                style(False),
            )

        has_clusters = ds.clusters is not None
        has_conditions = ds.conditions is not None
        has_samples = ds.samples is not None
        has_cell_types = ds.cell_types is not None

        embedding_flag = bool(getattr(profile, "embedding", False))
        dim_flag = embedding_flag
        colour_scale_flag = bool(getattr(profile, "colour_scale", False))

        options = []
        if getattr(profile, "split_by_condition", False):
            options.append({"label": "Split by condition", "value": OPT_SPLIT_BY_CONDITION})
        if getattr(profile, "is_3d", False) and embedding_flag:
            options.append({"label": "3D view", "value": OPT_IS_3D})
        show_options = bool(options)

        return (
            style(bool(getattr(profile, "clusters", True)) and has_clusters),
            style(bool(getattr(profile, "conditions", True)) and has_conditions),
            style(bool(getattr(profile, "samples", True)) and has_samples),
            style(bool(getattr(profile, "cell_types", True)) and has_cell_types),
            style(bool(getattr(profile, "genes", False))),
            style(embedding_flag),
            style(dim_flag),
            style(show_options),
            options,
            style(colour_scale_flag),
        )

    # ---------------------------------------------------------
    # Dimension selector population (labels from Dataset)
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.DIM_X, "options"),
        Output(IDs.Control.DIM_Y, "options"),
        Output(IDs.Control.DIM_Z, "options"),
        Output(IDs.Control.DIM_Z, "style"),
        Input(IDs.Control.EMBEDDING_SELECT, "value"),
        State(IDs.Control.DATASET_SELECT, "value"),
    )
    def update_dim_selectors(emb_key, dataset_name):
        empty: list[dict] = []
        hide_z = {"display": "none"}

        if not emb_key:
            return empty, empty, empty, hide_z

        ds = ctx.dataset_by_name.get(dataset_name) if dataset_name else None
        if ds is None:
            return empty, empty, empty, hide_z

        try:
            labels = ds.get_embedding_labels(emb_key)
        except KeyError:
            logger.warning(
                "Embedding key %r not found for dataset %r in update_dim_selectors",
                emb_key,
                dataset_name,
            )
            return empty, empty, empty, hide_z

        options = [{"label": labels[i], "value": i} for i in range(len(labels))]
        show_z = {} if len(labels) >= 3 else hide_z
        return options, options, options, show_z

    # ---------------------------------------------------------
    # Main figure: FilterState -> figure
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.MAIN_GRAPH, "figure"),
        Input(IDs.Store.FILTER_STATE, "data"),
    )
    def update_main_graph_from_state(fs_data):
        if fs_data is None:
            return _message_figure(
                "No dataset/view selected.",
                "Choose a dataset and view to see any plots.",
            )

        try:
            state = _parse_filter_state(fs_data)
        except Exception:
            logger.exception("Invalid filter state in main graph callback: %r", fs_data)
            return _error_figure("Internal error: invalid filter state.")

        ds = ctx.dataset_by_name.get(state.dataset_name)
        if ds is None:
            return _error_figure(
                f"The dataset '{state.dataset_name}' is not available. "
                "Try reloading the app or selecting a different dataset."
            )

        if fs_data.get("genes") and not state.genes:
            return _message_figure(
                "Selected genes not found in this dataset.",
                "Try a different gene symbol or clear the gene filter.",
            )

        if state.embedding and (state.dim_x is None or state.dim_y is None):
            return _message_figure(
                "Select dimensions for the embedding.",
                "Choose X and Y dimensions (and Z for 3D, if enabled) from the sidebar.",
            )

        try:
            view = ctx.registry.create(state.view_id, ds)

            logger.info(
                "update_main_graph",
                extra={
                    "view_id": state.view_id,
                    "dataset": state.dataset_name,
                    "n_clusters": len(state.clusters),
                    "n_conditions": len(state.conditions),
                    "n_samples": len(state.samples),
                    "n_cell_types": len(state.cell_types),
                    "n_genes": len(state.genes),
                },
            )

            data = view.timed_compute(state)

            if isinstance(data, pd.DataFrame) and data.empty:
                return _message_figure(
                    "No data to display.",
                    "Your current filters removed all cells for this view. "
                    "Try clearing one or more filters (clusters, conditions, samples, or cell types).",
                )

            if data is None:
                return _message_figure(
                    "No data returned by this view.",
                    "This can happen if there are no matching cells or the view "
                    "does not support the current settings.",
                )

            return view.render_figure(data, state)

        except Exception:
            logger.exception(
                "Error in update_main_graph_from_state",
                extra={"filter_state": fs_data},
            )
            return _error_figure(
                "The app hit an unexpected error. "
                "If this keeps happening, grab the logs and open an issue."
            )

    # ---------------------------------------------------------
    # Download CSV of current view data
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.DOWNLOAD_DATA, "data"),
        Input(IDs.Control.DOWNLOAD_DATA_BTN, "n_clicks"),
        State(IDs.Store.FILTER_STATE, "data"),
        prevent_initial_call=True,
    )
    def download_current_data(n_clicks, fs_data):
        if not n_clicks or fs_data is None:
            raise exceptions.PreventUpdate

        try:
            state = _parse_filter_state(fs_data)
        except Exception:
            raise exceptions.PreventUpdate

        ds = ctx.dataset_by_name.get(state.dataset_name)
        if ds is None:
            raise exceptions.PreventUpdate

        view = ctx.registry.create(state.view_id, ds)
        data = view.timed_compute(state)

        if not isinstance(data, pd.DataFrame) or data.empty:
            raise exceptions.PreventUpdate

        filename = f"{state.view_id}_{state.dataset_name.replace(' ', '_')}.csv"
        return dash_dcc.send_data_frame(data.to_csv, filename, index=False)

    # ---------------------------------------------------------
    # Save current figure -> metadata (via SessionService)
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
        n_clicks,
        fs_data,
        figure_label,
        session_data,
        active_session_id,
        active_figure_id,
    ):
        if not n_clicks:
            raise exceptions.PreventUpdate

        if fs_data is None:
            return session_data, active_session_id, "Failed to save: no current state.", active_figure_id

        try:
            state = _parse_filter_state(fs_data)
        except Exception:
            return session_data, active_session_id, "Failed to save: invalid state.", active_figure_id

        ds = ctx.dataset_by_name.get(state.dataset_name)
        if ds is None:
            return session_data, active_session_id, "Failed to save: dataset not found.", active_figure_id

        if active_session_id is None:
            active_session_id = generate_session_id()

        session = ensure_session(
            session_from_dict(session_data),
            session_id=active_session_id,
            app_version="0.0.0-dev",
            datasets_config_hash="unknown",
        )

        ds_key = getattr(ds, "key", ds.name)

        # IMPORTANT: store/persist dict only
        filter_state_dict = fs_data if isinstance(fs_data, dict) else {}

        session, new_active_figure_id, is_overwrite = save_figure(
            session,
            active_figure_id=active_figure_id,
            dataset_key=ds_key,
            view_id=state.view_id,
            filter_state=filter_state_dict,
            label=figure_label,
        )

        active_figure_id = new_active_figure_id
        view_label = _view_label(ctx, state.view_id)

        label_clean = None if figure_label is None else (str(figure_label).strip() or None)
        if is_overwrite:
            status = (
                f"Updated “{label_clean}” ({view_label}, {ds.name})."
                if label_clean
                else f"Updated {view_label} for {ds.name}."
            )
        else:
            status = (
                f"Saved “{label_clean}” ({view_label}, {ds.name})."
                if label_clean
                else f"Saved {view_label} for {ds.name}."
            )

        return session_to_dict(session), active_session_id, status, active_figure_id

    # ---------------------------------------------------------
    # Status bar
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.STATUS_BAR, "children"),
        Input(IDs.Store.FILTER_STATE, "data"),
    )
    def update_status_bar(fs_data):
        if fs_data is None:
            return html.Span(
                [
                    html.Strong("Dataset: "),
                    "No dataset selected",
                    " \u2022 ",
                    html.Strong("View: "),
                    "No view",
                    " \u2022 ",
                    html.Strong("Embedding: "),
                    "n/a",
                    " \u2022 ",
                    html.Strong("Genes: "),
                    "None",
                ]
            )

        try:
            state = _parse_filter_state(fs_data)
        except Exception:
            return html.Span(
                [
                    html.Strong("Dataset: "),
                    "Error",
                    " \u2022 ",
                    html.Strong("View: "),
                    "Error",
                    " \u2022 ",
                    html.Strong("Embedding: "),
                    "Error",
                    " \u2022 ",
                    html.Strong("Genes: "),
                    "Error",
                ]
            )

        ds = ctx.dataset_by_name.get(state.dataset_name)
        dataset_label = ds.name if ds is not None else f"{state.dataset_name} (missing)"
        view_label = _view_label(ctx, state.view_id) if state.view_id else "No view"
        emb_label = state.embedding or "Default embedding"

        genes = state.genes or []
        if not genes:
            gene_label = "None"
        elif len(genes) <= 3:
            gene_label = ", ".join(genes)
        else:
            gene_label = ", ".join(genes[:3]) + f" (+{len(genes) - 3} more)"

        return html.Span(
            [
                html.Strong("Dataset: "),
                dataset_label,
                " \u2022 ",
                html.Strong("View: "),
                view_label,
                " \u2022 ",
                html.Strong("Embedding: "),
                emb_label,
                " \u2022 ",
                html.Strong("Genes: "),
                gene_label,
            ]
        )

    # ---------------------------------------------------------
    # UI -> FilterState (canonical)  **USES ds.valid_sets()**
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Store.FILTER_STATE, "data"),
        Input(IDs.Control.DATASET_SELECT, "value"),
        Input(IDs.Control.VIEW_SELECT, "value"),
        Input(IDs.Control.CLUSTER_SELECT, "value"),
        Input(IDs.Control.CONDITION_SELECT, "value"),
        Input(IDs.Control.SAMPLE_SELECT, "value"),
        Input(IDs.Control.CELLTYPE_SELECT, "value"),
        Input(IDs.Control.GENE_SELECT, "value"),
        Input(IDs.Control.EMBEDDING_SELECT, "value"),
        Input(IDs.Control.DIM_X, "value"),
        Input(IDs.Control.DIM_Y, "value"),
        Input(IDs.Control.DIM_Z, "value"),
        Input(IDs.Control.OPTIONS_CHECKLIST, "value"),
        Input(IDs.Control.COLOUR_SCALE_SELECT, "value"),
        Input(IDs.Control.SAVED_FIGURE_LOAD_BTN, "n_clicks"),
    )
    def sync_filter_state_from_ui(
        dataset_name,
        view_id,
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
        colour_scale,
        _load_clicks,
    ):
        if dash.ctx.triggered_id == IDs.Control.SAVED_FIGURE_LOAD_BTN:
            raise exceptions.PreventUpdate

        if not dataset_name or not view_id:
            return None

        ds = ctx.dataset_by_name.get(dataset_name)
        if ds is None:
            return None

        valid = ds.valid_sets()

        orig_clusters = list(clusters or [])
        orig_conditions = list(conditions or [])
        orig_samples = list(samples or [])
        orig_cell_types = list(cell_types or [])
        orig_genes = list(genes or [])

        clusters = [str(c) for c in (clusters or []) if str(c) in valid.clusters]
        conditions = [str(c) for c in (conditions or []) if str(c) in valid.conditions]
        samples = [str(s) for s in (samples or []) if str(s) in valid.samples]
        cell_types = [str(ct) for ct in (cell_types or []) if str(ct) in valid.cell_types]
        genes = [g for g in (genes or []) if g in valid.genes]

        if (
            len(clusters) < len(orig_clusters)
            or len(conditions) < len(orig_conditions)
            or len(samples) < len(orig_samples)
            or len(cell_types) < len(orig_cell_types)
            or len(genes) < len(orig_genes)
        ):
            logger.info(
                "Dropped invalid filters when syncing UI for dataset %r: "
                "clusters %r -> %r, conditions %r -> %r, samples %r -> %r, "
                "cell_types %r -> %r, genes %r -> %r",
                dataset_name,
                orig_clusters,
                clusters,
                orig_conditions,
                conditions,
                orig_samples,
                samples,
                orig_cell_types,
                cell_types,
                orig_genes,
                genes,
            )

        if embedding and embedding not in ds.adata.obsm:
            embedding = ds.embedding_key if ds.embedding_key in ds.adata.obsm else None

        opts = set(options or [])
        split_by_condition = OPT_SPLIT_BY_CONDITION in opts
        is_3d = OPT_IS_3D in opts

        if embedding:
            try:
                labels = ds.get_embedding_labels(embedding)
            except KeyError:
                labels = []

            if labels:
                if dim_x is None and len(labels) >= 1:
                    dim_x = 0
                if dim_y is None and len(labels) >= 2:
                    dim_y = 1

                if is_3d:
                    if dim_z is None and len(labels) >= 3:
                        dim_z = 2
                    if dim_z is not None and dim_z >= len(labels):
                        dim_z = 2 if len(labels) >= 3 else None
                else:
                    dim_z = None

        state = FilterState(
            dataset_name=dataset_name,
            view_id=view_id,
            genes=genes,
            clusters=clusters,
            conditions=conditions,
            samples=samples,
            cell_types=cell_types,
            embedding=embedding,
            split_by_condition=split_by_condition,
            is_3d=is_3d,
            dim_x=dim_x,
            dim_y=dim_y,
            dim_z=dim_z,
            color_scale=colour_scale,
        )
        return state.to_dict()

    # ---------------------------------------------------------
    # Autosave / restore user state
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Store.USER_STATE, "data"),
        Input(IDs.Store.FILTER_STATE, "data"),
    )
    def autosave_user_state(fs_data):
        return fs_data or {}

    @app.callback(
        Output(IDs.Store.FILTER_STATE, "data", allow_duplicate=True),
        Input(IDs.Store.USER_STATE, "modified_timestamp"),
        State(IDs.Store.USER_STATE, "data"),
        prevent_initial_call=True,
    )
    def restore_user_state(_ts, state):
        if not state:
            raise exceptions.PreventUpdate

        try:
            restored = _parse_filter_state(state)
        except Exception:
            raise exceptions.PreventUpdate

        if restored.dataset_name not in ctx.dataset_by_name:
            logger.info("Skipping restore of user-state: dataset %r no longer available", restored.dataset_name)
            raise exceptions.PreventUpdate

        valid_view_ids = {cls.id for cls in ctx.registry.all_classes()}
        if restored.view_id not in valid_view_ids:
            logger.info("Skipping restore of user-state: view_id %r no longer registered", restored.view_id)
            raise exceptions.PreventUpdate

        return restored.to_dict()

    # ---------------------------------------------------------
    # Toggle Download button
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.DOWNLOAD_DATA_BTN, "disabled"),
        Input(IDs.Store.FILTER_STATE, "data"),
    )
    def toggle_download_button(fs_data):
        if fs_data is None:
            return True
        try:
            state = _parse_filter_state(fs_data)
        except Exception:
            return True
        if state.dataset_name not in ctx.dataset_by_name:
            return True
        return False

    # ---------------------------------------------------------
    # Saved figure dropdown (single authoritative callback)
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

        if current_value in ids or current_value == NEW_FIGURE_VALUE:
            value = current_value
        elif active_figure_id in ids:
            value = active_figure_id
        else:
            value = NEW_FIGURE_VALUE

        return options, value

    # ---------------------------------------------------------
    # Load a saved figure: apply full snapshot to UI + filter-state
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
        Output(IDs.Control.DIM_Z, "value", allow_duplicate=True),
        Output(IDs.Control.OPTIONS_CHECKLIST, "value", allow_duplicate=True),
        Output(IDs.Control.COLOUR_SCALE_SELECT, "value", allow_duplicate=True),
        Output(IDs.Control.FIGURE_LABEL_INPUT, "value", allow_duplicate=True),
        Output(IDs.Store.FILTER_STATE, "data", allow_duplicate=True),
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
            return (
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                "",
                no_update,
                None,
            )

        session = session_from_dict(session_data)
        if session is None:
            raise exceptions.PreventUpdate

        meta = next((f for f in session.figures if f.id == figure_id), None)
        if meta is None:
            raise exceptions.PreventUpdate

        ds = ctx.dataset_by_key.get(meta.dataset_key)

        # fallback 1: maybe dataset_key is actually the dataset name
        if ds is None:
            ds = ctx.dataset_by_name.get(meta.dataset_key)

        # fallback 2: match by name in case keys changed
        if ds is None:
            ds = next((d for d in ctx.datasets if d.name == meta.dataset_key), None)

        if ds is None:
            logger.warning(
                "Cannot load figure %r: dataset_key=%r not found in current datasets (%d available).",
                figure_id,
                meta.dataset_key,
                len(ctx.datasets),
            )
            raise exceptions.PreventUpdate

        fs_dict = meta.filter_state or {}
        try:
            state = _parse_filter_state(fs_dict)
        except Exception:
            logger.exception("Invalid filter_state on figure %r, aborting load", figure_id)
            raise exceptions.PreventUpdate

        state = replace(state, dataset_name=ds.name, view_id=meta.view_id)

        if state.embedding and state.embedding not in ds.adata.obsm:
            default_emb = getattr(ds, "embedding_key", None)
            state = replace(
                state,
                embedding=default_emb if (default_emb and default_emb in ds.adata.obsm) else None,
            )

        if state.embedding:
            try:
                labels = ds.get_embedding_labels(state.embedding)
            except KeyError:
                labels = []

            if labels:
                dx = state.dim_x if (state.dim_x is not None and 0 <= state.dim_x < len(labels)) else 0
                dy = state.dim_y if (state.dim_y is not None and 0 <= state.dim_y < len(labels)) else (
                    1 if len(labels) > 1 else 0
                )
                if state.is_3d:
                    dz = state.dim_z if (state.dim_z is not None and 0 <= state.dim_z < len(labels)) else (
                        2 if len(labels) > 2 else None
                    )
                else:
                    dz = None
                state = replace(state, dim_x=dx, dim_y=dy, dim_z=dz)

        options: list[str] = []
        if state.split_by_condition:
            options.append(OPT_SPLIT_BY_CONDITION)
        if state.is_3d:
            options.append(OPT_IS_3D)

        colour_scale = getattr(state, "color_scale", None)

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
            colour_scale,
            meta.label or "",
            state.to_dict(),
            figure_id,
        )
