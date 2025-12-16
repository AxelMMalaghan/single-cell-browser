from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import dash
from dash import Input, Output, State, html

from sc_browser.ui.helpers import get_filter_dropdown_options
from sc_browser.ui.ids import IDs

if TYPE_CHECKING:
    from sc_browser.ui.config import AppConfig

logger = logging.getLogger(__name__)

OPT_SPLIT_BY_CONDITION = "split_by_condition"
OPT_IS_3D = "is_3d"


def register_filter_callbacks(app: dash.Dash, ctx: AppConfig) -> None:
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
    # Live gene search (prefix, keep selected)
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
                style(True), style(True), style(True), style(True), style(True),
                style(True), style(True), style(True), options, style(False),
            )

        view = ctx.registry.create(view_id, ds)
        profile = getattr(view, "filter_profile", None)
        if profile is None:
            options = default_options()
            return (
                style(True), style(True), style(True), style(True), style(True),
                style(True), style(True), style(True), options, style(False),
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
    # Dimension selector population
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.DIM_X, "options"),
        Output(IDs.Control.DIM_Y, "options"),
        Output(IDs.Control.DIM_Z_SELECT, "options"),
        Output(IDs.Control.DIM_Z_CONTAINER, "style"),
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
            return empty, empty, empty, hide_z

        options = [{"label": labels[i], "value": i} for i in range(len(labels))]
        show_z = {} if len(labels) >= 3 else hide_z
        return options, options, options, show_z

