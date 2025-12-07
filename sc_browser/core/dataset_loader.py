from __future__ import annotations

from pathlib import Path

import anndata as ad

from sc_browser.config.model import DatasetConfig, ObsColumns
from sc_browser.core.dataset import Dataset


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

    # Semantic obs mapping from config
    obs_cols: ObsColumns = cfg.obs_columns

    group = cfg.raw.get("group", "Default")
    embedding_key = cfg.raw.get("embedding_key", "X_umap")

    # Pass cluster/condition explicitly and hand obs_cols to Dataset;
    # Dataset.__init__'s _normalise_obs_columns will accept ObsColumns.
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
