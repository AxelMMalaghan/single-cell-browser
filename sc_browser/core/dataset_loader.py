from __future__ import annotations
from pathlib import Path
from typing import Dict, Any

import anndata as ad

from sc_browser.config.inference_helpers import (
    _infer_cell_type_key,
    _infer_embedding_key,
    _infer_cluster_key,
    _infer_condition_key,
    _infer_sample_key,
)
from sc_browser.config.model import DatasetConfig
from sc_browser.core.dataset import Dataset


def from_config_entry(entry: Dict[str, Any]) -> Dataset:
    """
    Build a Dataset from a raw config entry (dict).

    This expects a mapping similar to what is found in JSON config files:
      - file/file_path/path: location of the .h5ad file
      - optional cluster_key / condition_key / embedding_key
      - optional obs_columns mapping (semantic name -> obs column)
    """
    file_value = entry.get("file_path") or entry.get("file") or entry.get("path")
    if file_value is None:
        raise KeyError("Config entry must contain 'file', 'file_path', or 'path'")

    file_path = Path(file_value)

    adata = ad.read_h5ad(file_path)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    # ---- obs column mapping ----
    obs_cols: Dict[str, str] = entry.get("obs_columns", {}) or {}

    # legacy top-level keys take precedence, then fall back to obs_columns
    cluster_key = entry.get("cluster_key") or obs_cols.get("cluster")
    condition_key = entry.get("condition_key") or obs_cols.get("condition")

    # Inference for missing cluster/condition
    if cluster_key is None:
        cluster_key = _infer_cluster_key(adata)
        if cluster_key is not None:
            obs_cols.setdefault("cluster", cluster_key)

    if condition_key is None:
        condition_key = _infer_condition_key(adata)
        if condition_key is not None:
            obs_cols.setdefault("condition", condition_key)

    # Optionally infer extra semantic columns
    sample_key = obs_cols.get("sample")
    if sample_key is None:
        inferred_sample = _infer_sample_key(adata)
        if inferred_sample is not None:
            obs_cols["sample"] = inferred_sample

    cell_type_key = obs_cols.get("cell_type")
    if cell_type_key is None:
        inferred_ct = _infer_cell_type_key(adata)
        if inferred_ct is not None:
            obs_cols["cell_type"] = inferred_ct

    # ---- embedding key ----
    embedding_key = entry.get("embedding_key")
    if embedding_key is None:
        embedding_key = _infer_embedding_key(adata)

    return Dataset(
        name=entry.get("name", "Unnamed dataset"),
        group=entry.get("group", "Default"),
        adata=adata,
        cluster_key=cluster_key,
        condition_key=condition_key,
        embedding_key=embedding_key,
        obs_columns=obs_cols,  # preserve inferred mapping
        file_path=file_path,
    )


def from_config(cfg: DatasetConfig) -> Dataset:
    """
    Build a Dataset from a DatasetConfig (new config layer).

    This is a thin adapter around `from_config_entry` so that all wiring
    logic lives in one place.
    """
    # Start from the raw JSON payload
    entry: Dict[str, Any] = dict(cfg.raw)

    # Normalise file path: our new schema uses "path", legacy might use "file"
    # from_config_entry expects one of "file", "file_path", or "path",
    # so we ensure "file_path" is present.
    entry.setdefault("file_path", str(cfg.path))

    # Ensure name/group are present, falling back to DatasetConfig
    entry.setdefault("name", cfg.name)
    entry.setdefault("group", entry.get("group", "Default"))

    return from_config_entry(entry)