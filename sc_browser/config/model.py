from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional


@dataclass(frozen=True)
class ObsColumns:
    """
    Semantic names for `.obs` columns used internally by the app.

    Purpose:
    - Provide a single, standardised object that describes "what column means what"
      for a given dataset.
    - Let the rest of the app talk in terms of semantic roles (cluster, condition,
      sample, etc.) instead of hardcoding raw `.obs` column names.

    Typical usage:
        ObsColumns(
            cluster="leiden",         # column in .obs holding cluster labels
            condition="treatment",    # column for experimental condition
            sample="sample_id",       # (optional) column for sample / donor id
        )

    Fields:
    - cluster: required; the column in `.obs` that encodes cluster / group labels
    - condition: optional; column encoding treatment / condition group
    - sample: optional; column encoding sample / donor / batch identifier
    """

    cluster: str
    condition: Optional[str] = None
    sample: Optional[str] = None


@dataclass
class DatasetConfig:
    """
    Parsed config entry for a single dataset.

    This is a thin wrapper around the raw config dict plus a few derived fields
    that make it easier to work with inside the app.

    Attributes
    ----------
    raw:
        The original config entry as loaded from JSON/YAML (one dataset block).
    source_path:
        Path to the config file this entry came from. Useful for error messages
        and relative paths.
    index:
        Index of this dataset entry within the config file (0-based), used to
        generate a fallback name if none is provided.
    """

    raw: Dict[str, Any]
    source_path: Path
    index: int

    @property
    def name(self) -> str:
        """
        Human-readable name for this dataset.

        Resolution order:
        - If `raw["name"]` exists, return that
        - Otherwise fall back to "Dataset <index>"
        """
        return self.raw.get("name", f"Dataset {self.index}")


@dataclass
class GlobalConfig:
    """
    Top-level config object for the application.

    This is the structured representation of your config file, after parsing.

    Attributes
    ----------
    ui_title:
        Title to display in the UI (e.g. app header).
    default_group:
        Default dataset "group" to select when the UI loads (e.g. a project or study).
    datasets:
        List of dataset configs available in this deployment.
    """

    ui_title: str
    default_group: str
    datasets: List[DatasetConfig]