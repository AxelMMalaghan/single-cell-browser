from __future__ import annotations

from typing import List, Optional, TYPE_CHECKING
import logging

import dash
import pandas as pd
import plotly.graph_objs as go
from dash import Input, Output, State

from sc_browser.core.filter_state import FilterState
from .helpers import get_filter_dropdown_options

from sc_browser.export.model import (
    FigureMetadata,
    new_session_metadata,
    session_from_dict,
    session_to_dict,
    generate_session_id,
    generate_figure_id,
    _now_iso,
)

if TYPE_CHECKING:
    from .context import AppContext   # type-only import

logger = logging.getLogger(__name__)


def _error_figure(message: str) -> go.Figure:
    """
    Generic fallback figure used when a view fails.
    """
    fig = go.Figure()
    fig.add_annotation(
        text=message,
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
        ds = ctx.dataset_by_name[dataset_name]
        name = ds.name
        meta = f"{ds.adata.n_obs} cells · {ds.adata.n_vars} genes"
        return name, meta

    # ---------------------------------------------------------
    # Reset filters when dataset changes
    # ---------------------------------------------------------
    @app.callback(
        Output("cluster-select", "value"),
        Output("condition-select", "value"),
        Output("sample-select", "value"),
        Output("celltype-select", "value"),
        Output("gene-select", "value"),
        Output("embedding-select", "value"),
        Output("dim-x-select", "value"),
        Output("dim-y-select", "value"),
        Output("dim-z-select", "value"),
        Input("dataset-select", "value"),
    )
    def reset_filters_on_dataset_change(dataset_name: str):
        ds = ctx.dataset_by_name[dataset_name]
        emb_val = ds.embedding_key if ds.embedding_key in ds.adata.obsm else None

        # Default dims: 0,1,(2)
        return [], [], [], [], [], emb_val, 0, 1, 2

    # ---------------------------------------------------------
    # Update dropdown options for filters
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
        ds = ctx.dataset_by_name[dataset_name]
        return get_filter_dropdown_options(ds)

    # ---------------------------------------------------------
    # Live gene search
    # ---------------------------------------------------------
    @app.callback(
        Output("gene-select", "options"),
        Input("gene-select", "search_value"),
        State("dataset-select", "value"),
        State("gene-select", "value"),
    )
    def update_gene_options(search_value, dataset_name, selected_genes):
        ds = ctx.dataset_by_name[dataset_name]
        all_genes = list(ds.genes)
        gene_set = set(all_genes)

        selected_genes = [g for g in (selected_genes or []) if g in gene_set]

        if not search_value or not search_value.strip():
            return [{"label": g, "value": g} for g in selected_genes]

        query = search_value.lower()
        matches = [g for g in all_genes if query in g.lower()]
        suggestions = matches[:100]

        union = list(dict.fromkeys(selected_genes + suggestions))
        return [{"label": g, "value": g} for g in union]

    # ---------------------------------------------------------
    # Hide/show filters based on active view
    # ---------------------------------------------------------
    # sc_browser/ui/callbacks_explore.py

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
        Input("view-tabs", "value"),
        State("dataset-select", "value"),
    )
    def update_filter_visibility(view_id: str, dataset_name: str):
        ds = ctx.dataset_by_name[dataset_name]
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

        # Build options list based on profile flags
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
    # Dimension selector population
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
        ds = ctx.dataset_by_name[dataset_name]

        coords = ds.get_embedding_matrix(emb_key)
        labels = ds.get_embedding_labels(emb_key)

        options = [{"label": labels[i], "value": i} for i in range(len(labels))]

        # Show Z selector only if 3 dims available
        show_z = {} if len(labels) >= 3 else {"display": "none"}

        return options, options, options, show_z

    # ---------------------------------------------------------
    # Main figure
    # ---------------------------------------------------------
    @app.callback(
        Output("main-graph", "figure"),
        Input("view-tabs", "value"),
        Input("dataset-select", "value"),
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
    def update_main_graph(
        view_id: str,
        dataset_name: str,
        clusters: List[str] | None,
        conditions: List[str] | None,
        samples: List[str] | None,
        cell_types: List[str] | None,
        genes: List[str] | None,
        embedding: Optional[str],
        dim_x: int,
        dim_y: int,
        dim_z: int,
        options: List[str] | None,
    ):
        try:
            ds = ctx.dataset_by_name[dataset_name]

            # Make sure genes exist in this dataset
            if genes:
                gene_set = set(ds.genes)
                genes = [g for g in genes if g in gene_set]

            state = FilterState(
                genes=genes or [],
                clusters=clusters or [],
                conditions=conditions or [],
                samples=samples or [],
                cell_types=cell_types or [],
                embedding=embedding,
                dim_x=dim_x,
                dim_y=dim_y,
                dim_z=dim_z,
                split_by_condition="split_by_condition" in (options or []),
                is_3d="is_3d" in (options or []),
            )

            view = ctx.registry.create(view_id, ds)

            logger.info(
                "update_main_graph",
                extra={
                    "view_id": view_id,
                    "dataset": dataset_name,
                    "n_clusters": len(state.clusters or []),
                    "n_conditions": len(state.conditions or []),
                    "n_genes": len(state.genes or []),
                },
            )

            data = view.timed_compute(state)
            fig = view.render_figure(data, state)
            return fig

        except Exception:
            # Log with context and return a friendly error figure
            logger.exception(
                "Error in update_main_graph",
                extra={
                    "view_id": view_id,
                    "dataset": dataset_name,
                    "clusters": clusters,
                    "conditions": conditions,
                    "samples": samples,
                    "cell_types": cell_types,
                    "genes": genes,
                    "embedding": embedding,
                    "dim_x": dim_x,
                    "dim_y": dim_y,
                    "dim_z": dim_z,
                    "options": options,
                },
            )
            return _error_figure("An error occurred while rendering this view. Check logs for details.")

    # ---------------------------------------------------------
    # Download CSV
    # ---------------------------------------------------------
    @app.callback(
        Output("download-data", "data"),
        Input("download-data-btn", "n_clicks"),
        State("view-tabs", "value"),
        State("dataset-select", "value"),
        State("cluster-select", "value"),
        State("condition-select", "value"),
        State("sample-select", "value"),
        State("celltype-select", "value"),
        State("gene-select", "value"),
        State("embedding-select", "value"),
        State("dim-x-select", "value"),
        State("dim-y-select", "value"),
        State("dim-z-select", "value"),
        State("options-checklist", "value"),
        prevent_initial_call=True,
    )
    def download_current_data(
        n_clicks,
        view_id: str,
        dataset_name: str,
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
        from dash import dcc as dash_dcc

        ds = ctx.dataset_by_name[dataset_name]
        state = FilterState(
            genes=genes or [],
            clusters=clusters or [],
            conditions=conditions or [],
            samples=samples or [],
            cell_types=cell_types or [],
            embedding=embedding,
            dim_x=dim_x,
            dim_y=dim_y,
            dim_z=dim_z,
            split_by_condition="split_by_condition" in (options or []),
        )

        view = ctx.registry.create(view_id, ds)
        data = view.timed_compute(state)
        if not isinstance(data, pd.DataFrame) or data.empty:
            return None

        filename = f"{view_id}_{dataset_name.replace(' ', '_')}.csv"
        return dash_dcc.send_data_frame(data.to_csv, filename, index=False)



    # ---------------------------------------------------------
    # Save current figure -> metadata + PNG
    # ---------------------------------------------------------
    @app.callback(
        Output("session-metadata", "data"),
        Output("active-session-id", "data"),
        Output("save-figure-status", "children"),
        Input("save-figure-btn", "n_clicks"),
        State("view-tabs", "value"),
        State("dataset-select", "value"),
        State("cluster-select", "value"),
        State("condition-select", "value"),
        State("sample-select", "value"),
        State("celltype-select", "value"),
        State("gene-select", "value"),
        State("embedding-select", "value"),
        State("dim-x-select", "value"),
        State("dim-y-select", "value"),
        State("dim-z-select", "value"),
        State("options-checklist", "value"),
        State("session-metadata", "data"),
        State("active-session-id", "data"),
        prevent_initial_call=True,
    )
    def save_current_figure(
        n_clicks,
        view_id: str,
        dataset_name: str,
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
        session_data,
        active_session_id,
    ):
        if not n_clicks:
            raise dash.exceptions.PreventUpdate

        # --- ensure we have a session id + session object ---
        if active_session_id is None:
            active_session_id = generate_session_id()

        session = session_from_dict(session_data)
        if session is None:
            # You can swap app_version / config_hash to anything real you have
            session = new_session_metadata(
                session_id=active_session_id,
                app_version="0.0.0-dev",
                datasets_config_hash="unknown",
            )

        # --- rebuild FilterState in exactly the same way as update_main_graph ---
        ds = ctx.dataset_by_name[dataset_name]

        if genes:
            gene_set = set(ds.genes)
            genes = [g for g in genes if g in gene_set]

        state = FilterState(
            genes=genes or [],
            clusters=clusters or [],
            conditions=conditions or [],
            samples=samples or [],
            cell_types=cell_types or [],
            embedding=embedding,
            dim_x=dim_x,
            dim_y=dim_y,
            dim_z=dim_z,
            split_by_condition="split_by_condition" in (options or []),
            is_3d="is_3d" in (options or []),
        )

        # dataset_key: for now use Dataset.name; later you can add a stable key field
        ds_key = getattr(ds, "key", ds.name)
        figure_id = generate_figure_id(session)

        meta = FigureMetadata.from_runtime(
            figure_id=figure_id,
            dataset_key=ds_key,
            view_id=view_id,
            state=state,
            view_params={},   # you can pass embedding/color_by here later if needed
            label=None,
            file_stem=None,
        )

        # --- render & write the image via export_service ---
        out_path = ctx.export_service.export_single(meta, session_id=session.session_id)
        meta.file_stem = out_path.stem

        # --- append to session + bump updated_at ---
        session.figures.append(meta)
        session.updated_at = _now_iso()

        status = f"Saved {meta.id} ({view_id}, {dataset_name})"
        return session_to_dict(session), active_session_id, status
