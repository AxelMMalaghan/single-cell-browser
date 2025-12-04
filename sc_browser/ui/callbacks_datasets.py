from __future__ import annotations

from typing import Optional, TYPE_CHECKING

import base64
import logging

import dash
from dash import Input, Output, State

from sc_browser.config.io import load_datasets
from sc_browser.config.write import save_dataset_config
from .helpers import obs_preview_table, dataset_status

if TYPE_CHECKING:
    from .context import AppContext   # type-only import


logger = logging.getLogger(__name__)


def register_dataset_callbacks(app: dash.Dash, ctx: "AppContext") -> None:
    @app.callback(
        Output("dm-status-text", "children"),
        Output("dm-summary-text", "children"),
        Output("dm-cluster-key", "options"),
        Output("dm-cluster-key", "value"),
        Output("dm-condition-key", "options"),
        Output("dm-condition-key", "value"),
        Output("dm-sample-key", "options"),
        Output("dm-sample-key", "value"),
        Output("dm-celltype-key", "options"),
        Output("dm-celltype-key", "value"),
        Output("dm-embedding-key", "options"),
        Output("dm-embedding-key", "value"),
        Output("dm-obs-preview", "children"),
        Output("dm-current-dataset", "children"),
        Input("dataset-select", "value"),
    )
    def update_dataset_manager(dataset_name: str):
        ds = ctx.dataset_by_name[dataset_name]

        status = dataset_status(ds)
        summary = f"{ds.adata.n_obs} cells · {ds.adata.n_vars} genes"

        obs_cols = list(ds.adata.obs.columns)
        obs_options = [{"label": c, "value": c} for c in obs_cols]

        cluster_val = ds.cluster_key if ds.cluster_key in ds.adata.obs.columns else None
        condition_val = ds.condition_key if ds.condition_key in ds.adata.obs.columns else None
        sample_val = ds.obs_columns.get("sample")
        celltype_val = ds.obs_columns.get("cell_type")

        emb_keys = list(ds.adata.obsm.keys())
        emb_options = [{"label": k, "value": k} for k in emb_keys]
        emb_val = ds.embedding_key if ds.embedding_key in ds.adata.obsm else (emb_keys[0] if emb_keys else None)

        obs_preview = obs_preview_table(ds)
        current_label = f"Current dataset: {ds.name}"

        return (
            status,
            summary,
            obs_options,
            cluster_val,
            obs_options,
            condition_val,
            obs_options,
            sample_val,
            obs_options,
            celltype_val,
            emb_options,
            emb_val,
            obs_preview,
            current_label,
        )

    @app.callback(
        Output("dm-save-status", "children"),
        Input("dm-save-btn", "n_clicks"),
        State("dataset-select", "value"),
        State("dm-cluster-key", "value"),
        State("dm-condition-key", "value"),
        State("dm-sample-key", "value"),
        State("dm-celltype-key", "value"),
        State("dm-embedding-key", "value"),
        prevent_initial_call=True,
    )
    def save_dataset_mapping(
        n_clicks: int,
        dataset_name: str,
        cluster_key: Optional[str],
        condition_key: Optional[str],
        sample_key: Optional[str],
        celltype_key: Optional[str],
        embedding_key: Optional[str],
    ):
        if not n_clicks:
            raise dash.exceptions.PreventUpdate

        ds = ctx.dataset_by_name[dataset_name]

        if cluster_key:
            ds.cluster_key = cluster_key
        if condition_key:
            ds.condition_key = condition_key
        if embedding_key:
            ds.embedding_key = embedding_key

        obs_cols = dict(ds.obs_columns)
        if sample_key:
            obs_cols["sample"] = sample_key
        if celltype_key:
            obs_cols["cell_type"] = celltype_key
        obs_cols["cluster"] = ds.cluster_key
        obs_cols["condition"] = ds.condition_key
        ds.obs_columns = obs_cols

        path = save_dataset_config(ds, ctx.config_root)
        return f"Saved mapping to {path.relative_to(ctx.config_root.parent)}"

    @app.callback(
        Output("dm-import-status", "children"),
        Output("dataset-select", "options"),
        Output("dataset-select", "value"),
        Input("dm-upload", "contents"),
        State("dm-upload", "filename"),
        State("dataset-select", "value"),
        prevent_initial_call=True,
    )
    def import_dataset(contents, filename, current_dataset_name):
        if not contents or not filename:
            raise dash.exceptions.PreventUpdate

        if not filename.endswith(".h5ad"):
            opts = [{"label": ds.name, "value": ds.name} for ds in ctx.datasets]
            current = current_dataset_name or (ctx.datasets[0].name if ctx.datasets else None)
            return (
                f"Unsupported file type for '{filename}'. Please upload a .h5ad file.",
                opts,
                current,
            )

        try:
            content_type, content_string = contents.split(",", 1)
            decoded = base64.b64decode(content_string)
        except Exception as e:
            logger.exception("Failed to decode uploaded file")
            opts = [{"label": ds.name, "value": ds.name} for ds in ctx.datasets]
            current = current_dataset_name or (ctx.datasets[0].name if ctx.datasets else None)
            return (f"Failed to decode uploaded file: {e}", opts, current)

        data_dir = ctx.config_root.parent / "data"
        data_dir.mkdir(parents=True, exist_ok=True)
        out_path = data_dir / filename

        try:
            with out_path.open("wb") as f:
                f.write(decoded)
        except Exception as e:
            logger.exception("Failed to write uploaded file to disk")
            opts = [{"label": ds.name, "value": ds.name} for ds in ctx.datasets]
            current = current_dataset_name or (ctx.datasets[0].name if ctx.datasets else None)
            return (f"Failed to save file to {out_path}: {e}", opts, current)

        try:
            ctx.global_config, new_datasets = load_datasets(ctx.config_root)
        except Exception as e:
            logger.exception("Reloading datasets after import failed")
            opts = [{"label": ds.name, "value": ds.name} for ds in ctx.datasets]
            current = current_dataset_name or (ctx.datasets[0].name if ctx.datasets else None)
            return (
                f"File saved to {out_path}, but reloading datasets failed: {e}",
                opts,
                current,
            )

        ctx.datasets = new_datasets
        ctx.dataset_by_name = {ds.name: ds for ds in ctx.datasets}

        imported_ds_name = None
        for ds in ctx.datasets:
            fp = getattr(ds, "file_path", None)
            if fp is not None and fp.name == filename:
                imported_ds_name = ds.name
                break

        if imported_ds_name is None and ctx.datasets:
            imported_ds_name = ctx.datasets[-1].name

        opts = [{"label": ds.name, "value": ds.name} for ds in ctx.datasets]
        selected = imported_ds_name or (ctx.datasets[0].name if ctx.datasets else None)

        status_msg = f"Imported '{filename}'"
        if imported_ds_name:
            status_msg += f" as dataset '{imported_ds_name}'."
        else:
            status_msg += " (could not match to a specific dataset; using existing list.)"

        return status_msg, opts, selected