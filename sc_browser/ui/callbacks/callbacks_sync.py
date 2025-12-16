from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

import dash
from dash import Input, Output, State, exceptions, html

from sc_browser.core.filter_state import FilterState
from sc_browser.ui.ids import IDs

if TYPE_CHECKING:
    from sc_browser.ui.config import AppConfig

logger = logging.getLogger(__name__)

OPT_SPLIT_BY_CONDITION = "split_by_condition"
OPT_IS_3D = "is_3d"


def _validate_and_build_state(ctx: AppConfig, inputs: dict[str, Any]) -> dict[str, Any] | None:
    """
    Pure helper to validate UI inputs against dataset capabilities and return a dict
    ready for FilterState.
    """
    # Extract values from args_grouping
    dataset_name = inputs.get("dataset", {}).get("value")
    view_id = inputs.get("view", {}).get("value")

    if not dataset_name or not view_id:
        return None

    ds = ctx.dataset_by_name.get(dataset_name)
    if ds is None:
        return None

    # Helper to safely get list values
    def get_list(key):
        val = inputs.get(key, {}).get("value")
        return list(val) if val else []

    # ---------------------------------------------------------
    # STABILITY FIX: Detect Trigger Source
    # ---------------------------------------------------------
    # If the dataset just changed, we MUST ignore the current values of dependent
    # filters (clusters, genes, etc.) because they are stale artifacts from the
    # previous dataset. The UI will clear them in a separate callback, but we
    # need the State to be clean *now* to prevent "ghost" states or name collisions.
    triggered_id = inputs.get("triggered_id")

    dataset_changed = triggered_id == IDs.Control.DATASET_SELECT

    if dataset_changed:
        # Force reset dependent filters
        clusters = []
        conditions = []
        samples = []
        cell_types = []
        genes = []
        # Embedding: attempt to keep if valid, or reset to default later
        embedding_val = inputs.get("embedding", {}).get("value")
    else:
        # Normal update: trust the inputs
        clusters = get_list("clusters")
        conditions = get_list("conditions")
        samples = get_list("samples")
        cell_types = get_list("cell_types")
        genes = get_list("genes")
        embedding_val = inputs.get("embedding", {}).get("value")

    # Validate against dataset (always safe)
    valid = ds.valid_sets()
    clusters = [str(c) for c in clusters if str(c) in valid.clusters]
    conditions = [str(c) for c in conditions if str(c) in valid.conditions]
    samples = [str(s) for s in samples if str(s) in valid.samples]
    cell_types = [str(ct) for ct in cell_types if str(ct) in valid.cell_types]
    genes = [g for g in genes if g in valid.genes]

    # Embedding logic
    # If dataset changed, we only keep embedding if it exists in new dataset,
    # otherwise default to the dataset's main embedding
    if dataset_changed:
        if embedding_val and embedding_val in ds.adata.obsm:
            embedding = embedding_val
        else:
            embedding = ds.embedding_key if ds.embedding_key in ds.adata.obsm else None
    else:
        # Standard validation
        embedding = embedding_val
        if embedding and embedding not in ds.adata.obsm:
            embedding = ds.embedding_key if ds.embedding_key in ds.adata.obsm else None

    # Options
    options = inputs.get("options", {}).get("value") or []
    opts_set = set(options)
    split_by_condition = OPT_SPLIT_BY_CONDITION in opts_set
    is_3d = OPT_IS_3D in opts_set

    # Dimensions
    # If dataset changed, reset dimensions to defaults (0,1,2) unless we want to get fancy
    if dataset_changed:
        dim_x, dim_y, dim_z = None, None, None
    else:
        dim_x = inputs.get("dim_x", {}).get("value")
        dim_y = inputs.get("dim_y", {}).get("value")
        dim_z = inputs.get("dim_z", {}).get("value")

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

    colour_scale = inputs.get("colour_scale", {}).get("value")

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


def register_sync_callbacks(app: dash.Dash, ctx: AppConfig) -> None:
    # ---------------------------------------------------------
    # UI -> FilterState (canonical)
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
        Input(IDs.Control.DIM_Z_SELECT, "value"),
        Input(IDs.Control.OPTIONS_CHECKLIST, "value"),
        Input(IDs.Control.COLOUR_SCALE_SELECT, "value"),
    )
    def sync_filter_state_from_ui(
            ds_val, view_val, clust_val, cond_val, samp_val, ct_val,
            gene_val, emb_val, dx_val, dy_val, dz_val, opt_val, col_val
    ):
        # Identify what triggered this callback
        # dash.ctx.triggered_id handles the ID resolution
        triggered_id = dash.ctx.triggered_id

        # We use ctx.args_grouping to pass named args to our helper
        inputs = {
            "triggered_id": triggered_id,  # Pass trigger info
            "dataset": {"value": ds_val},
            "view": {"value": view_val},
            "clusters": {"value": clust_val},
            "conditions": {"value": cond_val},
            "samples": {"value": samp_val},
            "cell_types": {"value": ct_val},
            "genes": {"value": gene_val},
            "embedding": {"value": emb_val},
            "dim_x": {"value": dx_val},
            "dim_y": {"value": dy_val},
            "dim_z": {"value": dz_val},
            "options": {"value": opt_val},
            "colour_scale": {"value": col_val},
        }

        return _validate_and_build_state(ctx, inputs)

    # ---------------------------------------------------------
    # Status Bar (Pure UI reflection of State)
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.STATUS_BAR, "children"),
        Input(IDs.Store.FILTER_STATE, "data"),
    )
    def update_status_bar(fs_data):
        if fs_data is None:
            return html.Span([html.Strong("Status: "), "No dataset selected"])

        try:
            # Quick parse, no strict object required for display
            ds_name = fs_data.get("dataset_name", "?")
            view_id = fs_data.get("view_id", "?")
            genes = fs_data.get("genes", [])
        except Exception:
            return html.Span([html.Strong("Status: "), "Error parsing state"])

        ds = ctx.dataset_by_name.get(ds_name)
        ds_label = ds.name if ds else ds_name

        # View label lookup
        view_label = view_id
        for cls in ctx.registry.all_classes():
            if cls.id == view_id:
                view_label = cls.label
                break

        gene_label = "None"
        if genes:
            gene_label = ", ".join(genes[:3]) + (f" (+{len(genes) - 3})" if len(genes) > 3 else "")

        return html.Span(
            [
                html.Strong("Dataset: "), ds_label, " \u2022 ",
                html.Strong("View: "), view_label, " \u2022 ",
                html.Strong("Genes: "), gene_label,
            ]
        )

    # ---------------------------------------------------------
    # Autosave / Restore User State
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
        if not state or not isinstance(state, dict):
            raise exceptions.PreventUpdate

        # Simple validation before restoring
        ds_name = state.get("dataset_name")
        if ds_name not in ctx.dataset_by_name:
            raise exceptions.PreventUpdate

        return state