from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad

from sc_browser.config.model import DatasetConfig, ObsColumns
from sc_browser.core.dataset import Dataset

logger = logging.getLogger(__name__)


class DatasetConfigError(ValueError):
    """
    Raised when a dataset config is structurally invalid for loading.
    """
    pass


def _ensure_unique_names(adata: ad.AnnData, cfg: DatasetConfig, path: Path) -> ad.AnnData:
    """
    Ensure obs_names and var_names are unique, logging what we do.
    """
    # obs_names
    if not adata.obs_names.is_unique:
        logger.warning(
            "Observation names are not unique for dataset '%s' (%s); "
            "calling .obs_names_make_unique() (in-memory fix)",
            cfg.name,
            path,
        )
        adata.obs_names_make_unique()

    # var_names
    if not adata.var_names.is_unique:
        logger.warning(
            "Variable names are not unique for dataset '%s' (%s); "
            "calling .var_names_make_unique() (in-memory fix)",
            cfg.name,
            path,
        )
        adata.var_names_make_unique()

    return adata


def _validate_obs_columns(adata: ad.AnnData, cfg: DatasetConfig, obs_cols: ObsColumns, path: Path) -> None:
    """
    Validate obs-column mappings **only if** they are explicitly configured.

    - cell_id: optional; if missing, we fall back to obs_names.
    - cluster / condition / sample / cell_type: optional; if provided but
      missing from .obs, that's a config error for curated datasets.
    """

    def check_optional(col_name: str | None, logical_name: str) -> None:
        if col_name is None:
            return
        if col_name not in adata.obs.columns:
            msg = (
                f"Dataset '{cfg.name}': obs_columns.{logical_name}='{col_name}' "
                f"not found in .obs"
            )
            logger.error(
                msg,
                extra={
                    "dataset": cfg.name,
                    "path": str(path),
                    logical_name: col_name,
                },
            )
            raise DatasetConfigError(msg)

    # cell_id is OPTIONAL – only complain if explicitly configured and wrong
    check_optional(obs_cols.cell_id, "cell_id")
    check_optional(obs_cols.cluster, "cluster")
    check_optional(obs_cols.condition, "condition")
    check_optional(obs_cols.sample, "sample")
    check_optional(obs_cols.cell_type, "cell_type")
    # batch is usually for integration / QC; also optional
    check_optional(obs_cols.batch, "batch")


def from_config(cfg: DatasetConfig) -> Dataset:
    """
    Materialise an AnnData-backed Dataset from a DatasetConfig.

    CRITICAL UPDATE:
    - Loads the .h5ad file in
    - This relies on the OS page cache rather than loading the full matrix into RAM.
    - Ensures obs_names and var_names are unique (in-memory).
    - Uses cfg.obs_columns (ObsColumns dataclass) to standardise obs column semantics.
    """
    path: Path = cfg.path
    if not path.is_file():
        raise DatasetConfigError(f"AnnData file not found at {path}")

    adata = ad.read_h5ad(path)
    adata = _ensure_unique_names(adata, cfg, path)

    # Semantic obs mapping from config; may be partially empty for discovered datasets
    obs_cols: ObsColumns = cfg.obs_columns

    # Only enforce columns that were *actually* specified in config
    _validate_obs_columns(adata, cfg, obs_cols, path)

    group = cfg.raw.get("group", "Default")
    embedding_key = cfg.raw.get("embedding_key", "X_umap")

    ds = Dataset(
        name=cfg.name,
        group=group,
        adata=adata,
        cluster_key=obs_cols.cluster,
        condition_key=obs_cols.condition,
        embedding_key=embedding_key,
        obs_columns=obs_cols,
        file_path=path,
    )

    logger.info(
        "Loaded dataset (loaded)",
        extra={
            "dataset": ds.name,
            "group": ds.group,
            "n_cells": ds.adata.n_obs,
            "n_genes": ds.adata.n_vars,
            "path": str(path),
        },
    )

    return ds

