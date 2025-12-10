from __future__ import annotations

from typing import Optional, TYPE_CHECKING

import base64
import logging

import dash
from dash import Input, Output, State

from sc_browser.config.loader import load_datasets
from sc_browser.config.new_config_writer import save_dataset_config
from .helpers import dataset_status

if TYPE_CHECKING:
    from .config import AppContext


from pathlib import Path

MAX_UPLOAD_BYTES = 1_000_000_000  # ~1 GB cap for .h5ad uploads; adjust as needed


logger = logging.getLogger(__name__)


def register_dataset_import_callbacks(app: dash.Dash, ctx: "AppContext") -> None:

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
        Output("dm-current-dataset", "children"),
        Input("dataset-select", "value"),
    )
    def update_dataset_mapping_and_status(dataset_name: str):
        # No dataset selected or not found – keep the UI stable and explain
        if not dataset_name or dataset_name not in ctx.dataset_by_name:
            empty_obs_options = []
            empty_emb_options = []
            return (
                "No dataset selected.",
                "",
                empty_obs_options,
                None,
                empty_obs_options,
                None,
                empty_obs_options,
                None,
                empty_obs_options,
                None,
                empty_emb_options,
                None,
                "Current dataset: none",
            )

        ds = ctx.dataset_by_name[dataset_name]

        status = dataset_status(ds)
        summary = f"{ds.adata.n_obs} cells · {ds.adata.n_vars} genes"

        obs_cols = list(ds.adata.obs.columns)
        obs_options = [{"label": c, "value": c} for c in obs_cols]

        cluster_val = ds.cluster_key if ds.cluster_key in obs_cols else None
        condition_val = ds.condition_key if ds.condition_key in obs_cols else None

        sample_val = ds.obs_columns.get("sample")
        if sample_val not in obs_cols:
            sample_val = None

        celltype_val = ds.obs_columns.get("cell_type")
        if celltype_val not in obs_cols:
            celltype_val = None

        emb_keys = list(ds.adata.obsm.keys())
        emb_options = [{"label": k, "value": k} for k in emb_keys]
        emb_val = (
            ds.embedding_key
            if ds.embedding_key in ds.adata.obsm
            else (emb_keys[0] if emb_keys else None)
        )

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

        if not dataset_name or dataset_name not in ctx.dataset_by_name:
            return "Cannot save mapping: no dataset selected."

        ds = ctx.dataset_by_name[dataset_name]

        # Track what we’re actually changing for a clearer message
        changed_fields = []

        if cluster_key and cluster_key != ds.cluster_key:
            ds.cluster_key = cluster_key
            changed_fields.append(f"cluster → {cluster_key}")
        if condition_key and condition_key != ds.condition_key:
            ds.condition_key = condition_key
            changed_fields.append(f"condition → {condition_key}")
        if embedding_key and embedding_key != ds.embedding_key:
            ds.embedding_key = embedding_key
            changed_fields.append(f"embedding → {embedding_key}")

        obs_cols = dict(ds.obs_columns)

        if sample_key and obs_cols.get("sample") != sample_key:
            obs_cols["sample"] = sample_key
            changed_fields.append(f"sample → {sample_key}")
        if celltype_key and obs_cols.get("cell_type") != celltype_key:
            obs_cols["cell_type"] = celltype_key
            changed_fields.append(f"cell_type → {celltype_key}")

        # Always keep these in sync with the dedicated keys
        obs_cols["cluster"] = ds.cluster_key
        obs_cols["condition"] = ds.condition_key
        ds.obs_columns = obs_cols

        # If literally nothing changed, don’t pretend we saved magic
        if not changed_fields:
            return f"No changes to save for dataset '{ds.name}'."

        path = save_dataset_config(ds, ctx.config_root)
        rel_path = path.relative_to(ctx.config_root.parent)

        changed_str = ", ".join(changed_fields)
        return (
            f"Updated mapping for '{ds.name}' "
            f"({changed_str}). Saved to {rel_path}."
        )

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

        # Base options in case we need to fall back
        def _current_options_and_value():
            opts = [{"label": ds.name, "value": ds.name} for ds in ctx.datasets]
            current = current_dataset_name or (ctx.datasets[0].name if ctx.datasets else None)
            return opts, current

        if not filename.endswith(".h5ad"):
            opts, current = _current_options_and_value()
            return (
                f"Unsupported file type for '{filename}'. "
                "Please upload a single-cell dataset in .h5ad format.",
                opts,
                current,
            )

        try:
            content_type, content_string = contents.split(",", 1)
            decoded = base64.b64decode(content_string)
        except Exception as e:
            logger.exception("Failed to decode uploaded file")
            opts, current = _current_options_and_value()
            return (
                f"Failed to read uploaded file '{filename}': {e}. "
                "Try re-uploading the .h5ad file.",
                opts,
                current,
            )

        # ---- SIZE GUARD ----
        if len(decoded) > MAX_UPLOAD_BYTES:
            logger.warning(
                "Rejecting uploaded file %r: size %d bytes exceeds MAX_UPLOAD_BYTES=%d",
                filename,
                len(decoded),
                MAX_UPLOAD_BYTES,
            )
            opts, current = _current_options_and_value()
            return (
                f"Upload too large for '{filename}' "
                f"({len(decoded) / (1024**2):.1f} MB; limit is {MAX_UPLOAD_BYTES / (1024**2):.0f} MB).",
                opts,
                current,
            )

        data_dir = ctx.config_root.parent / "data"
        data_dir.mkdir(parents=True, exist_ok=True)

        # ---- FILENAME SANITISATION ----
        # Strip any path components to avoid path traversal.
        safe_name = Path(filename).name
        if not safe_name.endswith(".h5ad"):
            safe_name = safe_name + ".h5ad"

        out_path = data_dir / safe_name

        try:
            with out_path.open("wb") as f:
                f.write(decoded)
        except Exception as e:
            logger.exception("Failed to write uploaded file to disk")
            opts, current = _current_options_and_value()
            return (
                f"Failed to save file to {out_path}: {e}. "
                "Check that the server has write permissions to the data directory.",
                opts,
                current,
            )

        try:
            ctx.global_config, new_datasets = load_datasets(ctx.config_root)
        except Exception as e:
            logger.exception("Reloading datasets after import failed")
            opts, current = _current_options_and_value()
            return (
                f"File '{safe_name}' was saved to {out_path}, "
                f"but reloading datasets failed: {e}. "
                "Fix the config or datasets and restart the app.",
                opts,
                current,
            )

        # Successful reload
        ctx.datasets = new_datasets
        ctx.dataset_by_name = {ds.name: ds for ds in ctx.datasets}

        imported_ds_name = None
        for ds in ctx.datasets:
            fp = getattr(ds, "file_path", None)
            if fp is not None and fp.name == safe_name:
                imported_ds_name = ds.name
                break

        # Fallback: if we can’t match by filename, guess “last dataset”
        if imported_ds_name is None and ctx.datasets:
            imported_ds_name = ctx.datasets[-1].name

        opts = [{"label": ds.name, "value": ds.name} for ds in ctx.datasets]
        selected = imported_ds_name or (ctx.datasets[0].name if ctx.datasets else None)

        if imported_ds_name:
            status_msg = (
                f"Imported '{safe_name}' as dataset '{imported_ds_name}'. "
                "It is now selected and can be configured in the Datasets → Import & mapping tab."
            )
        else:
            status_msg = (
                f"Imported '{safe_name}' and reloaded the dataset list, "
                "but could not confidently match it to a specific dataset entry. "
                "Check your config and dataset names."
            )

        return status_msg, opts, selected


