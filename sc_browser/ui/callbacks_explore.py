from __future__ import annotations

import logging

import dash
import pandas as pd
import plotly.graph_objs as go

from dash import Input, Output, State, exceptions, dcc as dash_dcc, html
from typing import Optional, TYPE_CHECKING


from .helpers import get_filter_dropdown_options
from sc_browser.core.filter_state import FilterState
from sc_browser.metadata_io.model import (
    FigureMetadata,
    new_session_metadata,
    session_from_dict,
    session_to_dict,
    generate_session_id,
    generate_figure_id,
    now_iso,
)

if TYPE_CHECKING:
    from .config import AppConfig

logger = logging.getLogger(__name__)


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
    return _message_figure(
        "Something went wrong while rendering this view.",
        details,
    )


def register_explore_callbacks(app: dash.Dash, ctx: "AppConfig") -> None:
    # ---------------------------------------------------------
    # Sidebar metadata
    # ---------------------------------------------------------
    @app.callback(
        Output("sidebar-dataset-name", "children"),
        Output("sidebar-dataset-meta", "children"),
        Input("dataset-select", "value"),
    )
    def update_sidebar_dataset_summary(dataset_name: str):
        ds = ctx.dataset_by_name.get(dataset_name)
        if ds is None:
            return "No dataset", "0 cells · 0 genes"
        return ds.name, f"{ds.adata.n_obs} cells · {ds.adata.n_vars} genes"

    # ---------------------------------------------------------
    # Update dropdown options for filters (per dataset)
    # ---------------------------------------------------------
    @app.callback(
        Output("cluster-select", "options"),
        Output("condition-select", "options"),
        Output("sample-select", "options"),
        Output("celltype-select", "options"),
        Output("embedding-select", "options"),
        Input("dataset-select", "value"),
    )
    def update_filters(dataset_name: str):
        ds = ctx.dataset_by_name.get(dataset_name)
        if ds is None:
            empty = []
            return empty, empty, empty, empty, empty
        return get_filter_dropdown_options(ds)

    # ---------------------------------------------------------
    # Live gene search (prefix, keep selected)
    # ---------------------------------------------------------
    @app.callback(
        Output("gene-select", "options"),
        Input("gene-select", "search_value"),
        Input("dataset-select", "value"),
        State("gene-select", "value"),
    )
    def update_gene_options(search_value, dataset_name, selected_genes):
        ds = ctx.dataset_by_name.get(dataset_name)
        if ds is None:
            return []

        all_genes = (
            list(getattr(ds, "genes", []))
            or list(getattr(ds, "gene_names", []))
            or list(ds.adata.var_names)
        )
        gene_set = set(all_genes)

        selected_genes = [g for g in (selected_genes or []) if g in gene_set]

        if not search_value or not str(search_value).strip():
            return [{"label": g, "value": g} for g in selected_genes]

        query = str(search_value).upper()
        matches = [g for g in all_genes if g.upper().startswith(query)]
        suggestions = matches[:50]

        union = list(dict.fromkeys(selected_genes + suggestions))
        return [{"label": g, "value": g} for g in union]

    # ---------------------------------------------------------
    # Hide/show filters based on active view
    # ---------------------------------------------------------
    @app.callback(
        Output("cluster-filter-container", "style"),
        Output("condition-filter-container", "style"),
        Output("sample-filter-container", "style"),
        Output("celltype-filter-container", "style"),
        Output("gene-filter-container", "style"),
        Output("embedding-filter-container", "style"),
        Output("dim-filter-container", "style"),
        Output("options-container", "style"),
        Output("options-checklist", "options"),
        Input("view-select", "value"),
        State("dataset-select", "value"),
    )
    def update_filter_visibility(view_id: str, dataset_name: str):
        ds = ctx.dataset_by_name.get(dataset_name)
        if ds is None:
            def style(_): return {}
            opts = [
                {"label": " Split by condition", "value": "split_by_condition"},
                {"label": " 3D view", "value": "is_3d"},
            ]
            return (
                style(True),
                style(True),
                style(True),
                style(True),
                style(True),
                style(True),
                style(True),
                style(True),
                opts,
            )

        view = ctx.registry.create(view_id, ds)
        profile = getattr(view, "filter_profile", None)

        def style(flag: bool) -> dict:
            return {} if flag else {"display": "none"}

        if profile is None:
            options = [
                {"label": " Split by condition", "value": "split_by_condition"},
                {"label": " 3D view", "value": "is_3d"},
            ]
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
            )

        embedding_flag = getattr(profile, "embedding", False)
        dim_flag = embedding_flag

        options = []
        if getattr(profile, "split_by_condition", False):
            options.append({"label": " Split by condition", "value": "split_by_condition"})
        if getattr(profile, "is_3d", False) and embedding_flag:
            options.append({"label": " 3D view", "value": "is_3d"})

        show_options = bool(options)

        return (
            style(getattr(profile, "clusters", True)),
            style(getattr(profile, "conditions", True)),
            style(getattr(profile, "samples", True)),
            style(getattr(profile, "cell_types", True)),
            style(getattr(profile, "genes", False)),
            style(embedding_flag),
            style(dim_flag),
            style(show_options),
            options,
        )

    # ---------------------------------------------------------
    # Dimension selector population (labels from Dataset)
    # ---------------------------------------------------------
    @app.callback(
        Output("dim-x-select", "options"),
        Output("dim-y-select", "options"),
        Output("dim-z-select", "options"),
        Output("dim-z-select", "style"),
        Input("embedding-select", "value"),
        State("dataset-select", "value"),
    )
    def update_dim_selectors(emb_key, dataset_name):
        empty = []
        hide_z = {"display": "none"}

        if not emb_key:
            return empty, empty, empty, hide_z

        ds = ctx.dataset_by_name.get(dataset_name)
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
        Output("main-graph", "figure"),
        Input("filter-state", "data"),
    )
    def update_main_graph_from_state(fs_data):
        if not fs_data:
            return _message_figure(
                "No dataset/view selected.",
                "Choose a dataset and view to see any plots.",
            )

        try:
            state = FilterState.from_dict(fs_data)
        except Exception:
            logger.exception("Invalid filter-state data: %r", fs_data)
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
        Output("download-data", "data"),
        Input("download-data-btn", "n_clicks"),
        State("filter-state", "data"),
        prevent_initial_call=True,
    )
    def download_current_data(n_clicks, fs_data):
        if not n_clicks or not fs_data:
            raise exceptions.PreventUpdate

        try:
            state = FilterState.from_dict(fs_data)
        except Exception:
            logger.exception("Invalid filter-state in download_current_data: %r", fs_data)
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
    # Save current figure -> metadata
    # ---------------------------------------------------------
    @app.callback(
        Output("session-metadata", "data"),
        Output("active-session-id", "data"),
        Output("save-figure-status", "children"),
        Output("active-figure-id", "data"),
        Input("save-figure-btn", "n_clicks"),
        State("filter-state", "data"),
        State("figure-label-input", "value"),
        State("session-metadata", "data"),
        State("active-session-id", "data"),
        State("active-figure-id", "data"),
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

        if not fs_data:
            return session_data, active_session_id, "Failed to save: no current state.", active_figure_id

        try:
            state = FilterState.from_dict(fs_data)
        except Exception:
            logger.exception("Invalid filter-state in save_current_figure: %r", fs_data)
            return session_data, active_session_id, "Failed to save: invalid state.", active_figure_id

        ds = ctx.dataset_by_name.get(state.dataset_name)
        if ds is None:
            return session_data, active_session_id, "Failed to save: dataset not found.", active_figure_id

        if active_session_id is None:
            active_session_id = generate_session_id()

        session = session_from_dict(session_data)
        if session is None:
            session = new_session_metadata(
                session_id=active_session_id,
                app_version="0.0.0-dev",
                datasets_config_hash="unknown",
            )

        if figure_label is None:
            label_clean = None
        else:
            stripped = figure_label.strip()
            label_clean = stripped or None

        ds_key = getattr(ds, "key", ds.name)

        existing_fig = None
        if active_figure_id:
            existing_fig = next((f for f in session.figures if f.id == active_figure_id), None)

        same_label_as_existing = (
            existing_fig is not None
            and (label_clean or "") == (existing_fig.label or "")
        )

        if existing_fig is not None and same_label_as_existing:
            figure_id = existing_fig.id
            is_overwrite = True
        else:
            figure_id = generate_figure_id(session)
            is_overwrite = False
            active_figure_id = figure_id

        meta = FigureMetadata.from_runtime(
            figure_id=figure_id,
            dataset_key=ds_key,
            view_id=state.view_id,
            filter_state=fs_data,
            view_params={},
            label=label_clean,
            file_stem=None,
        )

        if is_overwrite and existing_fig is not None:
            for idx, f in enumerate(session.figures):
                if f.id == existing_fig.id:
                    session.figures[idx] = meta
                    break
        else:
            session.figures.append(meta)

        session.updated_at = now_iso()

        view_label = _view_label(state.view_id)

        if is_overwrite:
            if label_clean:
                status = f"Updated “{label_clean}” ({view_label}, {ds.name})."
            else:
                status = f"Updated {view_label} for {ds.name}."
        else:
            if label_clean:
                status = f"Saved “{label_clean}” ({view_label}, {ds.name})."
            else:
                status = f"Saved {view_label} for {ds.name}."

        return session_to_dict(session), active_session_id, status, active_figure_id

    # ---------------------------------------------------------
    # Status bar
    # ---------------------------------------------------------
    @app.callback(
        Output("status-bar", "children"),
        Input("filter-state", "data"),
    )
    def update_status_bar(fs_data):
        if not fs_data:
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
            state = FilterState.from_dict(fs_data)
        except Exception:
            logger.exception("Invalid filter-state in update_status_bar: %r", fs_data)
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
        view_label = _view_label(state.view_id) if state.view_id else "No view"
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
    # UI -> FilterState (canonical)
    # ---------------------------------------------------------
    @app.callback(
        Output("filter-state", "data"),
        Input("dataset-select", "value"),
        Input("view-select", "value"),
        Input("cluster-select", "value"),
        Input("condition-select", "value"),
        Input("sample-select", "value"),
        Input("celltype-select", "value"),
        Input("gene-select", "value"),
        Input("embedding-select", "value"),
        Input("dim-x-select", "value"),
        Input("dim-y-select", "value"),
        Input("dim-z-select", "value"),
        Input("options-checklist", "value"),
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
    ):
        if not dataset_name or not view_id:
            return None

        ds = ctx.dataset_by_name.get(dataset_name)
        if ds is None:
            return None

        valid_clusters = set(ds.clusters.astype(str).unique())
        valid_conditions = (
            set(ds.conditions.astype(str).unique())
            if ds.conditions is not None
            else set()
        )
        valid_samples = (
            set(ds.samples.astype(str).unique())
            if getattr(ds, "samples", None) is not None
            else set()
        )
        valid_celltypes = (
            set(ds.cell_types.astype(str).unique())
            if getattr(ds, "cell_types", None) is not None
            else set()
        )
        valid_genes = (
            set(getattr(ds, "genes", []))
            or set(getattr(ds, "gene_names", []))
            or set(ds.adata.var_names)
        )

        orig_clusters = list(clusters or [])
        orig_conditions = list(conditions or [])
        orig_samples = list(samples or [])
        orig_cell_types = list(cell_types or [])
        orig_genes = list(genes or [])

        clusters = [str(c) for c in (clusters or []) if str(c) in valid_clusters]
        conditions = [str(c) for c in (conditions or []) if str(c) in valid_conditions]
        samples = [str(s) for s in (samples or []) if str(s) in valid_samples]
        cell_types = [str(ct) for ct in (cell_types or []) if str(ct) in valid_celltypes]
        genes = [g for g in (genes or []) if g in valid_genes]

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
        split_by_condition = "split_by_condition" in opts
        is_3d = "is_3d" in opts

        labels = []
        if embedding:
            try:
                labels = ds.get_embedding_labels(embedding)
            except KeyError:
                logger.warning(
                    "Embedding key %r not found for dataset %r while setting default dims",
                    embedding,
                    dataset_name,
                )

        if embedding and labels:
            if dim_x is None and len(labels) >= 1:
                dim_x = 0
            if dim_y is None and len(labels) >= 2:
                dim_y = 1
            if is_3d and dim_z is None and len(labels) >= 3:
                dim_z = 2
            if not is_3d:
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
        )

        return state.to_dict()

    # ---------------------------------------------------------
    # Autosave / restore user state
    # ---------------------------------------------------------
    @app.callback(
        Output("user-state", "data"),
        Input("filter-state", "data"),
    )
    def autosave_user_state(fs_data):
        return fs_data or {}

    @app.callback(
        Output("filter-state", "data", allow_duplicate=True),
        Input("user-state", "modified_timestamp"),
        State("user-state", "data"),
        prevent_initial_call=True,
    )
    def restore_user_state(_ts, state):
        if not state:
            raise exceptions.PreventUpdate

        try:
            restored = FilterState.from_dict(state)
        except Exception:
            logger.exception("Failed to restore user-state (invalid filter-state): %r", state)
            raise exceptions.PreventUpdate

        if restored.dataset_name not in ctx.dataset_by_name:
            logger.info(
                "Skipping restore of user-state: dataset %r no longer available",
                restored.dataset_name,
            )
            raise exceptions.PreventUpdate

        valid_view_ids = {cls.id for cls in ctx.registry.all_classes()}
        if restored.view_id not in valid_view_ids:
            logger.info(
                "Skipping restore of user-state: view_id %r no longer registered",
                restored.view_id,
            )
            raise exceptions.PreventUpdate

        return restored.to_dict()

    def _view_label(view_id: str) -> str:
        try:
            for cls in ctx.registry.all_classes():
                if getattr(cls, "id", None) == view_id:
                    return getattr(cls, "label", view_id)
        except Exception:
            pass
        return view_id

    # ---------------------------------------------------------
    # FilterState -> UI
    # ---------------------------------------------------------
    @app.callback(
        Output("dataset-select", "value", allow_duplicate=True),
        Output("view-select", "value", allow_duplicate=True),
        Output("cluster-select", "value", allow_duplicate=True),
        Output("condition-select", "value", allow_duplicate=True),
        Output("sample-select", "value", allow_duplicate=True),
        Output("celltype-select", "value", allow_duplicate=True),
        Output("gene-select", "value", allow_duplicate=True),
        Output("embedding-select", "value", allow_duplicate=True),
        Output("dim-x-select", "value", allow_duplicate=True),
        Output("dim-y-select", "value", allow_duplicate=True),
        Output("dim-z-select", "value", allow_duplicate=True),
        Input("filter-state", "data"),
        prevent_initial_call=True,
    )
    def sync_ui_from_filter_state(fs_data):
        if not fs_data:
            raise exceptions.PreventUpdate

        try:
            state = FilterState.from_dict(fs_data)
        except Exception:
            logger.exception(
                "Invalid filter-state in sync_ui_from_filter_state: %r",
                fs_data,
            )
            raise exceptions.PreventUpdate

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
        )

    # ---------------------------------------------------------
    # Toggle Download button
    # ---------------------------------------------------------
    @app.callback(
        Output("download-data-btn", "disabled"),
        Input("filter-state", "data"),
    )
    def toggle_download_button(fs_data):
        if not fs_data:
            return True

        try:
            state = FilterState.from_dict(fs_data)
        except Exception:
            return True

        if state.dataset_name not in ctx.dataset_by_name:
            return True

        return False

    # ---------------------------------------------------------
    # Saved figure dropdown (single authoritative callback)
    # ---------------------------------------------------------
    @app.callback(
        Output("saved-figure-select", "options"),
        Output("saved-figure-select", "value"),
        Input("session-metadata", "data"),
        State("saved-figure-select", "value"),
        State("active-figure-id", "data"),
    )
    def update_saved_figure_dropdown(session_data, current_value, active_figure_id):
        """
        Keep the saved-figure dropdown in sync with session metadata.

        - Rebuild options from session.figures
        - Try to keep the current selection if it's still valid
        - Otherwise, fall back to active_figure_id
        """
        session = session_from_dict(session_data)
        if session is None or not session.figures:
            return [], None

        opts = []
        ids = []
        for fig in session.figures:
            display = fig.label or fig.view_id
            display = f"{display} ({fig.dataset_key})"
            opts.append({"label": display, "value": fig.id})
            ids.append(fig.id)

        value = current_value if current_value in ids else None
        if value is None and active_figure_id in ids:
            value = active_figure_id

        return opts, value

    # ---------------------------------------------------------
    # Load a saved figure from the dropdown
    # ---------------------------------------------------------
    @app.callback(
        Output("dataset-select", "value", allow_duplicate=True),
        Output("active-figure-id", "data", allow_duplicate=True),
        Output("view-select", "value", allow_duplicate=True),
        Output("figure-label-input", "value", allow_duplicate=True),
        Output("cluster-select", "value", allow_duplicate=True),
        Output("condition-select", "value", allow_duplicate=True),
        Output("sample-select", "value", allow_duplicate=True),
        Output("celltype-select", "value", allow_duplicate=True),
        Output("gene-select", "value", allow_duplicate=True),
        Output("embedding-select", "value", allow_duplicate=True),
        Output("dim-x-select", "value", allow_duplicate=True),
        Output("dim-y-select", "value", allow_duplicate=True),
        Output("dim-z-select", "value", allow_duplicate=True),
        Output("options-checklist", "value", allow_duplicate=True),
        Input("saved-figure-select", "value"),
        State("session-metadata", "data"),
        prevent_initial_call=True,
    )
    def load_figure_from_dropdown(figure_id, session_data):
        if not figure_id or figure_id == "__new__":
            raise dash.exceptions.PreventUpdate

        session = session_from_dict(session_data)
        if session is None:
            raise dash.exceptions.PreventUpdate

        meta = next((f for f in session.figures if f.id == figure_id), None)
        if meta is None:
            raise dash.exceptions.PreventUpdate

        ds = ctx.dataset_by_key.get(meta.dataset_key)
        if ds is None:
            raise dash.exceptions.PreventUpdate

        fs = meta.filter_state or {}

        clusters = fs.get("clusters", [])
        conditions = fs.get("conditions", [])
        samples = fs.get("samples", [])
        cell_types = fs.get("cell_types", [])
        genes = fs.get("genes", [])

        embedding = fs.get("embedding", None)
        if not embedding or embedding not in ds.adata.obsm:
            default_emb = getattr(ds, "embedding_key", None)
            if default_emb and default_emb in ds.adata.obsm:
                embedding = default_emb
            else:
                embedding = None

        dim_x = fs.get("dim_x", 0)
        dim_y = fs.get("dim_y", 1)
        dim_z = fs.get("dim_z", 2)
        split_by_condition = fs.get("split_by_condition", False)
        is_3d = fs.get("is_3d", False)

        options = []
        if split_by_condition:
            options.append("split_by_condition")
        if is_3d:
            options.append("is_3d")

        dataset_name = ds.name

        return (
            dataset_name,
            figure_id,
            meta.view_id,
            meta.label or "",
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
        )
