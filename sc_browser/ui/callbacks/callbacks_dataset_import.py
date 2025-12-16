from __future__ import annotations

import base64
import logging
from pathlib import Path
from typing import Optional, TYPE_CHECKING

import dash
from dash import Input, Output, State

from sc_browser.config.dataset_loader import load_dataset_registry
from sc_browser.config.new_config_writer import save_dataset_config
from sc_browser.ui.helpers import dataset_status
from sc_browser.ui.ids import IDs
import dash_bootstrap_components as dbc
from dash import Input, Output, no_update

if TYPE_CHECKING:
    from sc_browser.ui.config import AppContext

MAX_UPLOAD_BYTES = 1_000_000_000
logger = logging.getLogger(__name__)


def _get_registry(ctx):
    """
    Best-effort access to the dataset registry in lazy mode.
    Adjust these attribute names once you know where you store it.
    """
    reg = getattr(ctx, "dataset_registry", None)
    if reg is not None:
        return reg

    for attr in ("registry", "dataset_registry", "datasets_registry"):
        reg = getattr(ctx.global_config, attr, None)
        if reg is not None:
            return reg

    return None


def _list_dataset_names(ctx) -> list[str]:
    datasets = getattr(ctx, "datasets", None)

    # Case 1: dict-like mapping (name -> Dataset/handle/config)
    if isinstance(datasets, dict):
        return list(datasets.keys())

    # Case 2: sequence (list/tuple)
    if isinstance(datasets, (list, tuple)):
        if not datasets:
            return []
        first = datasets[0]
        if isinstance(first, str):
            return list(datasets)
        if hasattr(first, "name"):
            return [ds.name for ds in datasets]

    # Case 3: registry object reachable via ctx/global_config
    reg = _get_registry(ctx)
    if reg is not None:
        # Prefer an explicit dataset_names list if it exists
        for attr in ("dataset_names", "names"):
            val = getattr(reg, attr, None)
            if callable(val):
                return list(val())
            if isinstance(val, (list, tuple)):
                return list(val)

        # Sometimes the registry stores configs in a dict
        for attr in ("datasets", "dataset_configs", "by_name", "configs"):
            val = getattr(reg, attr, None)
            if isinstance(val, dict):
                return list(val.keys())

    # Fallback: whatever is currently materialised
    by_name = getattr(ctx, "dataset_by_name", None)
    if isinstance(by_name, dict):
        return list(by_name.keys())

    return []


def _materialise_dataset(ctx, name: str):
    ds = ctx.dataset_by_name.get(name)
    if ds is not None:
        return ds

    reg = _get_registry(ctx)

    # If ctx.datasets is a dict and already holds Dataset objects, use it
    datasets = getattr(ctx, "datasets", None)
    if isinstance(datasets, dict):
        maybe = datasets.get(name)
        if maybe is not None and hasattr(maybe, "name") and hasattr(maybe, "adata"):
            ds = maybe
            ctx.dataset_by_name[name] = ds
            ctx.dataset_by_key[getattr(ds, "key", ds.name)] = ds
            if ctx.default_dataset is None:
                ctx.default_dataset = ds
            return ds

    if reg is None:
        return None

    for meth in ("materialise", "materialize", "load", "get", "get_dataset", "dataset"):
        fn = getattr(reg, meth, None)
        if callable(fn):
            try:
                ds = fn(name)
            except Exception:
                logger.exception("Dataset registry method %s(%r) failed", meth, name)
                ds = None

            if ds is not None:
                ctx.dataset_by_name[name] = ds
                ctx.dataset_by_key[getattr(ds, "key", ds.name)] = ds
                if ctx.default_dataset is None:
                    ctx.default_dataset = ds
                return ds

    return None



def _dataset_select_options(ctx) -> list[dict]:
    names = _list_dataset_names(ctx)
    return [{"label": n, "value": n} for n in names]


def _coerce_selected_name(ctx, current_name: str | None) -> str | None:
    names = _list_dataset_names(ctx)
    if not names:
        return None
    if current_name and current_name in names:
        return current_name
    return names[0]


def register_dataset_import_callbacks(app: dash.Dash, ctx: "AppContext") -> None:
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
            ds = ctx.dataset_by_name.get(dataset_name)
            if ds is None:
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
            ds = ctx.dataset_by_name.get(dataset_name)
            if ds is None:
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
    # Import dataset upload (.h5ad), reload registry, refresh dropdown
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
                f"Unsupported file type for '{filename}'. "
                "Please upload a single-cell dataset in .h5ad format.",
                opts,
                current,
            )

        try:
            _content_type, content_string = contents.split(",", 1)
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
            ctx.global_config, reg_or_ds = load_dataset_registry(ctx.config_root)
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

        # reg_or_ds in lazy mode is NOT necessarily a list.
        # It might be a dict or a registry-like object.
        ctx.datasets = reg_or_ds

        # Clear materialised caches; they may now be stale
        ctx.dataset_by_name = {}
        ctx.dataset_by_key = {}
        ctx.default_dataset = None

        names = _list_dataset_names(ctx)
        opts = [{"label": n, "value": n} for n in names]

        # Best-effort: select imported dataset by filename stem in dataset name
        stem = Path(safe_name).stem.lower()
        imported_ds_name = next((n for n in names if stem in n.lower()), None)

        selected = imported_ds_name or (
            current_dataset_name if current_dataset_name in names else (names[-1] if names else None)
        )

        if imported_ds_name:
            status_msg = (
                f"Imported '{safe_name}'. "
                f"Dataset list reloaded and '{imported_ds_name}' is now selected."
            )
        else:
            status_msg = (
                f"Imported '{safe_name}' and reloaded the dataset list, "
                "but could not confidently match it to a specific dataset entry. "
                "Select it from the dataset dropdown."
            )

        return status_msg, opts, selected


    @app.callback(
        Output(IDs.Control.DATASET_TOAST, "is_open"),
        Output(IDs.Control.DATASET_TOAST, "children"),
        Output(IDs.Control.DATASET_TOAST, "icon"),
        Output(IDs.Control.DATASET_STATUS_BADGE, "children"),
        Input(IDs.Control.DATASET_SELECT, "value"),
        prevent_initial_call=True,
    )
    def dataset_ui_feedback(dataset_name: str | None):
        if not dataset_name:
            badge = dbc.Badge("No dataset selected", color="secondary", className="ms-2")
            return True, "No dataset selected.", "secondary", badge

        # If you have lazy mode, you can show whether it’s materialised
        is_materialised = dataset_name in getattr(ctx, "dataset_by_name", {})

        if is_materialised:
            msg = f"'{dataset_name}' is loaded."
            icon = "success"
            badge = dbc.Badge(f"Loaded: {dataset_name}", color="success", className="ms-2")
        else:
            msg = f"Selected '{dataset_name}'. (Lazy mode: will load when needed.)"
            icon = "info"
            badge = dbc.Badge(f"Selected: {dataset_name}", color="info", className="ms-2")

        return True, msg, icon, badge