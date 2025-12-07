from __future__ import annotations

from pathlib import Path

import anndata as ad
import logging

from sc_browser.config.model import DatasetConfig, ObsColumns
from sc_browser.core.dataset import Dataset

logger = logging.getLogger(__name__)


def from_config(cfg: DatasetConfig) -> Dataset:
    """
    Materialise an AnnData-backed Dataset from a DatasetConfig.

    - Loads the .h5ad file from cfg.path
    - Uses cfg.obs_columns (ObsColumns dataclass) to standardise obs column semantics
    - Applies per-dataset settings (group, embedding key) from cfg.raw
    """
    path: Path = cfg.path
    if not path.is_file():
        raise FileNotFoundError(f"AnnData file not found at {path}")

    adata = ad.read_h5ad(path)

    if not adata.obs_names.is_unique:
        logger.warning(""
                       "Observation names are not unique for dataset '%s' (%s), "
                       "calling .obs_names_make_unique()",
                       cfg.name,
                       path)
        adata = adata.copy()
        adata.obs_names_make_unique()


    # Semantic obs mapping from config
    obs_cols: ObsColumns = cfg.obs_columns

    group = cfg.raw.get("group", "Default")
    embedding_key = cfg.raw.get("embedding_key", "X_umap")

    dataset = Dataset(
        name=cfg.name,
        group=group,
        adata=adata,
        cluster_key=obs_cols.cluster,
        condition_key=obs_cols.condition,
        embedding_key=embedding_key,
        obs_columns=obs_cols,
        file_path=path,)


    logger.info(
        "Loaded dataset",
        extra={
            "dataset": dataset.name,
            "group": dataset.group,
            "n_cells": dataset.adata.n_obs,
            "n_genes": dataset.adata.n_vars,
            "path": str(path),
        },
    )

    return dataset

