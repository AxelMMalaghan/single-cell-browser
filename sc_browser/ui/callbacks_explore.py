from __future__ import annotations

from typing import List, Optional, TYPE_CHECKING
import logging

import dash
import pandas as pd
import plotly.graph_objs as go
from dash import Input, Output, State

from sc_browser.core.filter_state import FilterState
from .helpers import get_filter_dropdown_options

if TYPE_CHECKING:
    from .context import AppContext   # type-only import

logger = logging.getLogger(__name__)


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
    @app.callback(
        Output("cluster-filter-container", "style"),
        Output("condition-filter-container", "style"),
        Output("sample-filter-container", "style"),
        Output("celltype-filter-container", "style"),
        Output("gene-filter-container", "style"),
        Output("embedding-filter-container", "style"),
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
            return (style(True),) * 6

        return (
            style(getattr(profile, "clusters", True)),
            style(getattr(profile, "conditions", True)),
            style(getattr(profile, "samples", True)),
            style(getattr(profile, "cell_types", True)),
            style(getattr(profile, "genes", False)),
            style(getattr(profile, "embedding", False)),
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

        view = ctx.registry.create(view_id, ds)

        try:
            data = view.compute_data(state)
            fig = view.render_figure(data, state)
            return fig

        except Exception as e:
            print("DEBUG ERROR:", repr(e))
            raise

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
        data = view.compute_data(state)
        if not isinstance(data, pd.DataFrame) or data.empty:
            return None

        filename = f"{view_id}_{dataset_name.replace(' ', '_')}.csv"
        return dash_dcc.send_data_frame(data.to_csv, filename, index=False)