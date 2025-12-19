from __future__ import annotations

import logging
import os
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
    if not adata.obs_names.is_unique:
        logger.warning(
            "Observation names are not unique for dataset '%s' (%s); "
            "calling .obs_names_make_unique() (in-memory fix)",
            cfg.name,
            path,
        )
        adata.obs_names_make_unique()

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
    Validate obs-column mappings only if they are explicitly configured.
    """

    def check_optional(col_name: str | None, logical_name: str) -> None:
        if col_name is None:
            return
        if col_name not in adata.obs.columns:
            msg = (
                f"Dataset '{cfg.name}': obs_columns.{logical_name}='{col_name}' "
                f"not found in .obs"
            )
            logger.error(msg, extra={"dataset": cfg.name, "path": str(path), logical_name: col_name})
            raise DatasetConfigError(msg)

    check_optional(obs_cols.cell_id, "cell_id")
    check_optional(obs_cols.cluster, "cluster")
    check_optional(obs_cols.condition, "condition")
    check_optional(obs_cols.sample, "sample")
    check_optional(obs_cols.cell_type, "cell_type")
    check_optional(obs_cols.batch, "batch")


def from_config(cfg: DatasetConfig) -> Dataset:
    """
    Materialise an AnnData-backed Dataset from a DatasetConfig.
    """
    path: Path = cfg.path

    # RESOLUTION LOGIC: Use SC_BROWSER_DATA_ROOT for relative paths
    if not path.is_absolute():
        data_root = os.environ.get("SC_BROWSER_DATA_ROOT")
        if data_root:
            root_path = Path(data_root)
            resolved_path = root_path / path

            # Fallback for redundant 'data/' prefix
            if not resolved_path.is_file() and path.parts[0] == "data":
                alt_path = root_path / Path(*path.parts[1:])
                if alt_path.is_file():
                    resolved_path = alt_path
            path = resolved_path

    if not path.is_file():
        raise DatasetConfigError(f"AnnData file not found at {path}.")

    adata = ad.read_h5ad(path)
    adata = _ensure_unique_names(adata, cfg, path)

    obs_cols: ObsColumns = cfg.obs_columns
    _validate_obs_columns(adata, cfg, obs_cols, path)

    group = cfg.raw.get("group", "Default")
    embedding_key = cfg.raw.get("embedding_key", "X_umap")

    return Dataset(
        name=cfg.name,
        group=group,
        adata=adata,
        cluster_key=obs_cols.cluster,
        condition_key=obs_cols.condition,
        embedding_key=embedding_key,
        obs_columns=obs_cols,
        file_path=path,
    )
