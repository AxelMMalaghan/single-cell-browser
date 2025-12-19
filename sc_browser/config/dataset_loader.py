from __future__ import annotations

import logging
import os
from pathlib import Path

import anndata as ad

from sc_browser.core.configs import DatasetConfig
from sc_browser.core.dataset import Dataset

logger = logging.getLogger(__name__)


class DatasetConfigError(ValueError):
    """Raised when a dataset config is structurally invalid for loading."""
    pass


def _ensure_unique_names(adata: ad.AnnData, cfg: DatasetConfig, path: Path) -> ad.AnnData:
    if not adata.obs_names.is_unique:
        logger.warning(f"Making obs_names unique for {cfg.name}")
        adata.obs_names_make_unique()
    if not adata.var_names.is_unique:
        logger.warning(f"Making var_names unique for {cfg.name}")
        adata.var_names_make_unique()
    return adata


def from_config(cfg: DatasetConfig) -> Dataset:
    """Materialise an AnnData-backed Dataset from a DatasetConfig."""
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
        raise DatasetConfigError(f"AnnData file not found at {path}")

    logger.info(f"Loading AnnData from {path}")
    adata = ad.read_h5ad(path)
    adata = _ensure_unique_names(adata, cfg, path)

    group = cfg.raw.get("group", "Default")
    embedding_key = cfg.raw.get("embedding_key", "X_umap")

    return Dataset(
        name=cfg.name,
        group=group,
        adata=adata,
        cluster_key=cfg.obs_columns.cluster,
        condition_key=cfg.obs_columns.condition,
        embedding_key=embedding_key,
        obs_columns=cfg.obs_columns,
        file_path=path,
    )