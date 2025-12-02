from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

import anndata
import pandas as pd

# ---- Updated Semantic Mapping and DatasetConfig ----



@dataclass
class DatasetConfig:
    """
    Parsed config entry for a single dataset.
    """
    raw: Dict[str, Any]
    source_path: Path
    index: int

    @property
    def name(self) -> str:
        return self.raw.get("name", f"Dataset {self.index}")

    @property
    def path(self) -> Path:
        return Path(self.raw["file"])

    @property
    def obs_columns(self) -> ObsColumns:
        return ObsColumns(**self.raw.get("obs_columns", {}))

    @classmethod
    def from_raw(cls, raw: Dict[str, Any], source_path: Path, index: int) -> DatasetConfig:
        return cls(raw=raw, source_path=source_path, index=index)


@dataclass
class GlobalConfig:
    ui_title: str
    default_group: str
    datasets: List[DatasetConfig]


@dataclass(frozen=True)
class ObsColumns:
    """
    Semantic names for `.obs` columns used internally by the app.
    """
    cell_id: str
    cluster: Optional[str] = None
    condition: Optional[str] = None
    sample: Optional[str] = None
    batch: Optional[str] = None
    cell_type: Optional[str] = None


# ---- Preview + Inference Functions ----

def infer_column_roles(obs_df: pd.DataFrame) -> Dict[str, Optional[str]]:
    """Heuristic-based inference of semantic column roles."""
    colnames = obs_df.columns
    roles = {
        "cluster": next((col for col in colnames if col.lower() in ["cluster", "leiden", "louvain"]), None),
        "condition": next((col for col in colnames if col.lower() in ["condition", "treatment", "group"]), None),
        "sample": next((col for col in colnames if col.lower() in ["sample", "donor", "replicate"]), None),
        "batch": next((col for col in colnames if col.lower() in ["batch", "batch_id"]), None),
        "cell_type": next((col for col in colnames if col.lower() in ["cell_type", "type", "annotation"]), None),
    }
    return roles

def preview_adata_columns(file_path: Path) -> Dict[str, Any]:
    adata = anndata.read_h5ad(file_path, backed='r')
    obs_df = adata.obs
    inferred_roles = infer_column_roles(obs_df)
    return {
        "obs_columns": list(obs_df.columns),
        "suggested_roles": inferred_roles,
        "n_obs": adata.n_obs,
        "n_vars": adata.n_vars,
        "obsm_keys": list(adata.obsm.keys()),
    }

