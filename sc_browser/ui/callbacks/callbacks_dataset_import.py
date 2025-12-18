from __future__ import annotations

import base64
import binascii
import errno
import logging
from pathlib import Path
from typing import Optional, TYPE_CHECKING, List

import dash
from dash import Input, Output, State

from sc_browser.config.dataset_loader import load_dataset_registry
from sc_browser.config.new_config_writer import save_dataset_config
from sc_browser.config.save_ds_config import write_dataset_config_for_uploaded_file
from sc_browser.ui.helpers import dataset_status
from sc_browser.ui.ids import IDs
from sc_browser.services.dataset_service import DatasetManager

from sc_browser.core.inference import (
    infer_cluster_key,
    infer_condition_key,
    infer_sample_key,
    infer_cell_type_key,
    infer_embedding_key
)

if TYPE_CHECKING:
    from sc_browser.ui.config import AppConfig

MAX_UPLOAD_BYTES = 1_000_000_000
logger = logging.getLogger(__name__)


def _list_dataset_names(ctx: AppConfig) -> List[str]:
    return sorted(list(ctx.dataset_by_name.keys()))


def _materialise_dataset(ctx: AppConfig, name: str):
    return ctx.dataset_by_name.get(name)


def _dataset_select_options(ctx: AppConfig) -> list[dict]:
    names = _list_dataset_names(ctx)
    return [{"label": n, "value": n} for n in names]


def _coerce_selected_name(ctx: AppConfig, current_name: str | None) -> str | None:
    names = _list_dataset_names(ctx)
    if not names:
        return None
    if current_name and current_name in names:
        return current_name
    return names[0]


def register_dataset_import_callbacks(app: dash.Dash, ctx: AppConfig) -> None:
    # ---------------------------------------------------------
    # Mapping/status panel updates when dataset-select changes
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.DM_STATUS_TEXT, "children"),
        Output(IDs.Control.DM_SUMMARY_TEXT, "children"),
        Output(IDs.Control.DM_CLUSTER_KEY, "options"),
        Output(IDs.Control.DM_CLUSTER_KEY, "value"),
        Output(IDs.Control.DM_CONDITION_KEY, "options"),
        Output(IDs.Control.DM_CONDITION_KEY, "value"),
        Output(IDs.Control.DM_SAMPLE_KEY, "options"),
        Output(IDs.Control.DM_SAMPLE_KEY, "value"),
        Output(IDs.Control.DM_CELLTYPE_KEY, "options"),
        Output(IDs.Control.DM_CELLTYPE_KEY, "value"),
        Output(IDs.Control.DM_EMBEDDING_KEY, "options"),
        Output(IDs.Control.DM_EMBEDDING_KEY, "value"),
        Output(IDs.Control.DM_CURRENT_DATASET, "children"),
        Input(IDs.Control.DATASET_SELECT, "value"),
    )
    def update_dataset_mapping_and_status(dataset_name: str | None):
        ds = None
        if dataset_name:
            ds = _materialise_dataset(ctx, dataset_name)

        if ds is None:
            empty_obs_options: list[dict] = []
            empty_emb_options: list[dict] = []
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

        status = dataset_status(ds)
        summary = f"{ds.adata.n_obs} cells · {ds.adata.n_vars} genes"

        obs_cols = list(ds.adata.obs.columns)
        obs_options = [{"label": c, "value": c} for c in obs_cols]

        # --- CLUSTER ---
        cluster_val = ds.cluster_key if ds.cluster_key in obs_cols else None
        if cluster_val is None:
            cluster_val = infer_cluster_key(ds.adata)

        # --- CONDITION ---
        condition_val = ds.condition_key if ds.condition_key in obs_cols else None
        if condition_val is None:
            condition_val = infer_condition_key(ds.adata)

        # --- SAMPLE ---
        sample_val = ds.obs_columns.get("sample")
        if sample_val not in obs_cols:
            sample_val = infer_sample_key(ds.adata)

        # --- CELL TYPE ---
        celltype_val = ds.obs_columns.get("cell_type")
        if celltype_val not in obs_cols:
            celltype_val = infer_cell_type_key(ds.adata)

        # --- EMBEDDING ---
        emb_keys = list(ds.adata.obsm.keys())
        emb_options = [{"label": k, "value": k} for k in emb_keys]

        emb_val = ds.embedding_key
        if emb_val not in ds.adata.obsm:
            emb_val = infer_embedding_key(ds.adata)

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

    # ---------------------------------------------------------
    # Save mapping changes
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.DM_SAVE_STATUS, "children"),
        Input(IDs.Control.DM_SAVE_BTN, "n_clicks"),
        State(IDs.Control.DATASET_SELECT, "value"),
        State(IDs.Control.DM_CLUSTER_KEY, "value"),
        State(IDs.Control.DM_CONDITION_KEY, "value"),
        State(IDs.Control.DM_SAMPLE_KEY, "value"),
        State(IDs.Control.DM_CELLTYPE_KEY, "value"),
        State(IDs.Control.DM_EMBEDDING_KEY, "value"),
        prevent_initial_call=True,
    )
    def save_dataset_mapping(
            n_clicks: int,
            dataset_name: str | None,
            cluster_key: Optional[str],
            condition_key: Optional[str],
            sample_key: Optional[str],
            celltype_key: Optional[str],
            embedding_key: Optional[str],
    ):
        if not n_clicks:
            raise dash.exceptions.PreventUpdate

        ds = None
        if dataset_name:
            ds = _materialise_dataset(ctx, dataset_name)

        if ds is None:
            return "Cannot save mapping: no dataset selected."

        changed_fields: list[str] = []

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

        # Keep semantic keys aligned
        obs_cols["cluster"] = ds.cluster_key
        obs_cols["condition"] = ds.condition_key
        ds.obs_columns = obs_cols

        refresh = getattr(ds, "refresh_obs_indexes", None)
        if callable(refresh):
            refresh()
        else:
            logger.warning(
                "Dataset.refresh_obs_indexes() not found on %r. "
                "Mapping changes may not be reflected until restart.",
                ds.name,
            )

        if not changed_fields:
            return f"No changes to save for dataset '{ds.name}'."

        path = save_dataset_config(ds, ctx.config_root)
        rel_path = path.relative_to(ctx.config_root.parent)

        changed_str = ", ".join(changed_fields)

        return f"Updated mapping for '{ds.name}' ({changed_str}). Saved to {rel_path}."

    # ---------------------------------------------------------
    # Import dataset upload
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.DM_IMPORT_STATUS, "children"),
        Output(IDs.Control.DATASET_SELECT, "options"),
        Output(IDs.Control.DATASET_SELECT, "value"),
        Input(IDs.Control.DM_UPLOAD, "contents"),
        State(IDs.Control.DM_UPLOAD, "filename"),
        State(IDs.Control.DATASET_SELECT, "value"),
        prevent_initial_call=True,
    )
    def import_dataset(contents, filename, current_dataset_name):
        if not contents or not filename:
            raise dash.exceptions.PreventUpdate

        def _current_options_and_value():
            opts = _dataset_select_options(ctx)
            current = _coerce_selected_name(ctx, current_dataset_name)
            return opts, current

        if not filename.endswith(".h5ad"):
            opts, current = _current_options_and_value()
            return (
                f"Unsupported file type for '{filename}'. Please upload a .h5ad file.",
                opts,
                current,
            )

        # 1. DECODING
        try:
            _content_type, content_string = contents.split(",", 1)
            decoded = base64.b64decode(content_string)
        except (ValueError, binascii.Error) as e:
            logger.error("Corrupted upload data for %s: %s", filename, e)
            opts, current = _current_options_and_value()
            return f"The uploaded file '{filename}' appears to be corrupted.", opts, current

        if len(decoded) > MAX_UPLOAD_BYTES:
            opts, current = _current_options_and_value()
            return f"File '{filename}' exceeds the 1GB limit.", opts, current

        # 2. PREPARE PATHS
        data_dir = ctx.config_root.parent / "data"
        data_dir.mkdir(parents=True, exist_ok=True)
        safe_name = Path(filename).name
        if not safe_name.endswith(".h5ad"): safe_name += ".h5ad"
        out_path = data_dir / safe_name

        # 3. ATOMIC WRITE
        try:
            with out_path.open("xb") as f:
                f.write(decoded)
        except FileExistsError:
            opts, current = _current_options_and_value()
            return f"A dataset named '{safe_name}' already exists on the server.", opts, current
        except PermissionError:
            logger.error("Permission denied writing to %s", out_path)
            opts, current = _current_options_and_value()
            return "Server Error: Permission denied when saving the file.", opts, current
        except OSError as e:
            opts, current = _current_options_and_value()
            if e.errno == errno.ENOSPC:
                return "Server Error: The disk is full. Cannot save dataset.", opts, current
            logger.exception("Disk I/O error during import")
            return f"Server Error: Failed to save file ({e.strerror})", opts, current

        # 4. REGISTRATION & RELOAD
        try:
            write_dataset_config_for_uploaded_file(
                config_root=ctx.config_root,
                dataset_name=Path(safe_name).stem,
                file_path=out_path,
                group="Imported"
            )
            ctx.global_config, reg_or_ds = load_dataset_registry(ctx.config_root)
        except Exception as e:
            logger.exception("Post-upload registration failed")
            opts, current = _current_options_and_value()
            return f"File saved, but failed to register dataset: {e}", opts, current

        # Refresh Manager
        if isinstance(ctx.dataset_by_name, DatasetManager):
            ctx.dataset_by_name.refresh_config(reg_or_ds)

        names = _list_dataset_names(ctx)
        opts = [{"label": n, "value": n} for n in names]
        stem = Path(safe_name).stem.lower()
        imported_ds_name = next((n for n in names if stem in n.lower()), None)
        selected = imported_ds_name or (
            current_dataset_name if current_dataset_name in names else (names[-1] if names else None))

        return f"Successfully imported '{safe_name}'.", opts, selected