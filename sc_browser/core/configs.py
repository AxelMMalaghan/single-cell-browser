from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional


@dataclass(frozen=True)
class ObsColumns:
    """
    Semantic names for `.obs` columns used internally by the app.

    All fields are optional so that auto-discovered datasets without explicit
    obs_columns config still parse. Views must be defensive when fields are None.
    """

    cell_id: Optional[str] = None
    cluster: Optional[str] = None
    condition: Optional[str] = None
    sample: Optional[str] = None
    batch: Optional[str] = None
    cell_type: Optional[str] = None


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
        """
        Return the .h5ad path for this dataset.

        Supports both:
        - new schema:  "path": "data/foo.h5ad"
        - legacy:      "file": "data/foo.h5ad" or "file_path": "data/foo.h5ad"
        """
        raw_path = (
            self.raw.get("path") or self.raw.get("file") or self.raw.get("file_path")
        )
        if raw_path is None:
            raise KeyError(
                f"No 'path', 'file', or 'file_path' in dataset config: {self.raw}"
            )
        return Path(raw_path)

    @property
    def obs_columns(self) -> ObsColumns:
        """
        Return semantic obs column mapping as an ObsColumns instance.

        If 'obs_columns' is missing or partial, unspecified fields default to None.
        """
        raw_cols = self.raw.get("obs_columns", {})
        return ObsColumns(**raw_cols)

    @classmethod
    def from_raw(
        cls, raw: Dict[str, Any], source_path: Path, index: int
    ) -> DatasetConfig:
        return cls(raw=raw, source_path=source_path, index=index)


@dataclass
class GlobalConfig:
    ui_title: str
    default_group: str
    datasets: List[DatasetConfig]
    data_root: Optional[Path] = None
