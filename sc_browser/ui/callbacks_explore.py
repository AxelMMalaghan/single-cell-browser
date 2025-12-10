from __future__ import annotations

import logging
from typing import Optional, TYPE_CHECKING

import dash
import pandas as pd
import plotly.graph_objs as go
from dash import Input, Output, State, exceptions, dcc as dash_dcc, html

from .helpers import get_filter_dropdown_options



from sc_browser.core.filter_state import FilterState
from sc_browser.metadata_io.model import (
    FigureMetadata,
    new_session_metadata,
    session_from_dict,
    session_to_dict,
    generate_session_id,
    generate_figure_id,
    _now_iso,
)

if TYPE_CHECKING:
    from .config import AppContext

logger = logging.getLogger(__name__)


def _message_figure(title: str, details: Optional[str] = None) -> go.Figure:
    """
    Generic figure used to show user-facing status / error messages.
    """
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
    """
    Convenience wrapper for error situations.
    """
    return _message_figure(
        "Something went wrong while rendering this view.",
        details,
    )


def register_explore_callbacks(app: dash.Dash, ctx: "AppContext") -> None:
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

        # Use the same helper as layout/_build_filter_panel
        return get_filter_dropdown_options(ds)

    # ---------------------------------------------------------
    # Live gene search (prefix, keep selected)
    # ---------------------------------------------------------
    @app.callback(
        Output("gene-select", "options"),
        Input("gene-select", "search_value"),
        State("dataset-select", "value"),
        State("gene-select", "value"),
    )
    def update_gene_options(search_value, dataset_name, selected_genes):
        ds = ctx.dataset_by_name.get(dataset_name)
        if ds is None:
            raise exceptions.PreventUpdate

        # dataset API: prefer ds.genes, fall back to ds.gene_names or var_names
        all_genes = (
            list(getattr(ds, "genes", []))
            or list(getattr(ds, "gene_names", []))
            or list(ds.adata.var_names)
        )
        gene_set = set(all_genes)

        selected_genes = [g for g in (selected_genes or []) if g in gene_set]

        if not search_value or not search_value.strip():
            # Only show selected when there's no query
            return [{"label": g, "value": g} for g in selected_genes]

        query = search_value.upper()
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
            # no dataset → just show everything
            style = lambda _: {}
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
            # Fallback: show everything, both toggles
            options = [
                {"label": " Split by condition", "value": "split_by_condition"},
                {"label": " 3D view", "value": "is_3d"},
            ]
            return (
                style(True),  # clusters
                style(True),  # conditions
                style(True),  # samples
                style(True),  # cell_types
                style(True),  # genes
                style(True),  # embedding
                style(True),  # dims
                style(True),  # options-container
                options,
            )

        embedding_flag = getattr(profile, "embedding", False)
        dim_flag = embedding_flag  # dims only when embedding is active

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
        if not emb_key:
            empty = []
            hide_z = {"display": "none"}
            return empty, empty, empty, hide_z

        ds = ctx.dataset_by_name.get(dataset_name)
        if ds is None:
            empty = []
            hide_z = {"display": "none"}
            return empty, empty, empty, hide_z

        try:
            labels = ds.get_embedding_labels(emb_key)
        except KeyError:
            empty = []
            hide_z = {"display": "none"}
            logger.warning(
                "Embedding key %r not found for dataset %r in update_dim_selectors",
                emb_key,
                dataset_name,
            )
            return empty, empty, empty, hide_z

        options = [{"label": labels[i], "value": i} for i in range(len(labels))]
        show_z = {} if len(labels) >= 3 else {"display": "none"}
        return options, options, options, show_z

    # ---------------------------------------------------------
    # Main figure
    # ---------------------------------------------------------

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
            # User had gene selections, but none are valid in this dataset
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

            fig = view.render_figure(data, state)
            return fig

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
            # Nothing sensible to download
            raise exceptions.PreventUpdate

        filename = f"{state.view_id}_{state.dataset_name.replace(' ', '_')}.csv"
        return dash_dcc.send_data_frame(data.to_csv, filename, index=False)

    # ---------------------------------------------------------
    # Save current figure -> metadata + PNG
    # ---------------------------------------------------------

    @app.callback(
        Output("session-metadata", "data"),
        Output("active-session-id", "data"),
        Output("save-figure-status", "children"),
        Input("save-figure-btn", "n_clicks"),
        State("filter-state", "data"),
        State("figure-label-input", "value"),
        State("session-metadata", "data"),
        State("active-session-id", "data"),
        prevent_initial_call=True,
    )
    def save_current_figure(
        n_clicks,
        fs_data,
        figure_label,
        session_data,
        active_session_id,
    ):
        if not n_clicks:
            raise exceptions.PreventUpdate

        if not fs_data:
            return session_data, active_session_id, "Failed to save: no current state."

        try:
            state = FilterState.from_dict(fs_data)
        except Exception:
            logger.exception("Invalid filter-state in save_current_figure: %r", fs_data)
            return session_data, active_session_id, "Failed to save: invalid state."

        ds = ctx.dataset_by_name.get(state.dataset_name)
        if ds is None:
            return session_data, active_session_id, "Failed to save: dataset not found."

        if active_session_id is None:
            active_session_id = generate_session_id()

        session = session_from_dict(session_data)
        if session is None:
            session = new_session_metadata(
                session_id=active_session_id,
                app_version="0.0.0-dev",
                datasets_config_hash="unknown",
            )

        # clean / normalise label
        if figure_label is None:
            label_clean = None
        else:
            stripped = figure_label.strip()
            label_clean = stripped or None

        ds_key = getattr(ds, "key", ds.name)
        figure_id = generate_figure_id(session)

        # IMPORTANT: write filter_state as a dict, consistent with reports code
        meta = FigureMetadata.from_runtime(
            figure_id=figure_id,
            dataset_key=ds_key,
            view_id=state.view_id,
            filter_state=state.to_dict(),
            view_params={},
            label=label_clean,
            file_stem=None,
        )

        out_path = ctx.export_service.export_single(meta, session_id=session.session_id)
        meta.file_stem = out_path.stem

        session.figures.append(meta)
        session.updated_at = _now_iso()

        view_label = _view_label(state.view_id)

        if label_clean:
            status = f"Saved “{label_clean}” ({view_label}, {ds.name})."
        else:
            status = f"Saved {view_label} for {ds.name}."

        return session_to_dict(session), active_session_id, status

    # ---------------------------------------------------------
    # Status bar
    # ---------------------------------------------------------

    @app.callback(
        Output("status-bar", "children"),
        Input("filter-state", "data"),
    )
    def update_status_bar(fs_data):
        """
        Status bar derived from the canonical FilterState.

        This ensures the embedding / genes we show actually match the
        normalised state (after validation, dataset switching, etc).
        """
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
        """
        Build the canonical FilterState from current UI selections.

        This is the ONLY place that writes filter-state.data.
        """
        if not dataset_name or not view_id:
            # No meaningful state yet
            return None

        ds = ctx.dataset_by_name.get(dataset_name)
        if ds is None:
            return None

        # ---- validate + normalise selections against dataset ----
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

        # Keep originals for logging
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

        # Log if we dropped anything
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

        # Embedding sanity: if invalid, prefer dataset default or None
        if embedding and embedding not in ds.adata.obsm:
            embedding = ds.embedding_key if ds.embedding_key in ds.adata.obsm else None

        opts = set(options or [])
        split_by_condition = "split_by_condition" in opts
        is_3d = "is_3d" in opts

        # Let dims be None until explicitly chosen, so the main-graph guard works
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
    # Autosave user state (FilterState-based)
    # ---------------------------------------------------------
    @app.callback(
        Output("user-state", "data"),
        Input("filter-state", "data"),
    )
    def autosave_user_state(fs_data):
        """
        Persist the last known FilterState for restoration.
        """
        return fs_data or {}

    # ---------------------------------------------------------
    # Restore user state on reload -> FilterState
    # ---------------------------------------------------------

    @app.callback(
        Output("filter-state", "data", allow_duplicate=True),
        Input("user-state", "modified_timestamp"),
        State("user-state", "data"),
        prevent_initial_call=True,
    )
    def restore_user_state(_ts, state):
        """
        Restore the previous FilterState when the app reloads.

        We validate it against current datasets/views to avoid resurrecting
        completely invalid state after config changes.
        """
        if not state:
            raise exceptions.PreventUpdate

        try:
            restored = FilterState.from_dict(state)
        except Exception:
            logger.exception("Failed to restore user-state (invalid filter-state): %r", state)
            raise exceptions.PreventUpdate

        # dataset must still exist
        if restored.dataset_name not in ctx.dataset_by_name:
            logger.info(
                "Skipping restore of user-state: dataset %r no longer available",
                restored.dataset_name,
            )
            raise exceptions.PreventUpdate

        # view_id must still be registered
        valid_view_ids = {cls.id for cls in ctx.registry.all_classes()}
        if restored.view_id not in valid_view_ids:
            logger.info(
                "Skipping restore of user-state: view_id %r no longer registered",
                restored.view_id,
            )
            raise exceptions.PreventUpdate

        return restored.to_dict()


    def _view_label(view_id: str) -> str:
        """
        Resolve a human-readable view label from the registry.
        Falls back to the raw id if not found.
        """
        try:
            for cls in ctx.registry.all_classes():
                if getattr(cls, "id", None) == view_id:
                    return getattr(cls, "label", view_id)
        except Exception:
            pass
        return view_id

    @app.callback(
        Output("dataset-select", "value", allow_duplicate=True),
        Output("view-select", "value", allow_duplicate=True),
        Output("cluster-select", "value", allow_duplicate=True),
        Output("condition-select", "value", allow_duplicate=True),
        Output("sample-select", "value", allow_duplicate=True),
        Output("celltype-select", "value", allow_duplicate=True),
        Output("gene-select", "value", allow_duplicate=True),
        Output("embedding-select", "value", allow_duplicate=True),
        Input("filter-state", "data"),
        prevent_initial_call=True,
    )
    def sync_ui_from_filter_state(fs_data):
        if not fs_data:
            raise exceptions.PreventUpdate

        try:
            state = FilterState.from_dict(fs_data)
        except Exception:
            logger.exception("Invalid filter-state in sync_ui_from_filter_state: %r", fs_data)
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
        )


    @app.callback(
        Output("download-data-btn", "disabled"),
        Input("filter-state", "data"),
    )
    def toggle_download_button(fs_data):
        """
        Disable the Download CSV button when there is no valid filter state.
        """
        if not fs_data:
            return True

        try:
            state = FilterState.from_dict(fs_data)
        except Exception:
            return True

        if state.dataset_name not in ctx.dataset_by_name:
            return True

        # For now, assume all views *may* return tabular data.
        # If you later add a flag like view.supports_table_export, check it here.
        return False




