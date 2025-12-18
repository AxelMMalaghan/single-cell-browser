from __future__ import annotations

import uuid
import logging
from dataclasses import dataclass, field
from typing import Any, Optional, Dict, List
from datetime import datetime, timezone

# Re-import FilterState for proper reconstruction
from sc_browser.core.filter_state import FilterState

logger = logging.getLogger(__name__)


# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

def generate_figure_id() -> str:
    """Generate a globally unique, immutable ID for figures."""
    return f"fig-{uuid.uuid4().hex[:12]}"


def generate_session_id() -> str:
    """Generate a globally unique, immutable ID for sessions."""
    return f"session-{uuid.uuid4().hex[:8]}"


def now_iso() -> str:
    """Return a current UTC timestamp in ISO-8601 format."""
    return datetime.now(timezone.utc).isoformat()


# -------------------------------------------------------------------------
# Per-figure metadata
# -------------------------------------------------------------------------

@dataclass
class FigureMetadata:
    """
    Description of a single saved figure.
    Ensures filter_state is a FilterState object, not just a dict.
    """
    id: str
    dataset_key: str
    view_id: str
    filter_state: FilterState
    view_params: Dict[str, Any] = field(default_factory=dict)
    label: Optional[str] = None
    file_stem: Optional[str] = None
    created_at: str = field(default_factory=now_iso)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> FigureMetadata:
        """Reconstructs FigureMetadata, ensuring FilterState is an object."""
        fs_raw = data.get("filter_state") or {}
        # FIX: Reconstitute FilterState if it is a dictionary
        fs_obj = fs_raw if isinstance(fs_raw, FilterState) else FilterState.from_dict(fs_raw)

        return cls(
            id=data.get("id", generate_figure_id()),
            dataset_key=data["dataset_key"],
            view_id=data["view_id"],
            filter_state=fs_obj,
            view_params=data.get("view_params", {}),
            label=data.get("label"),
            file_stem=data.get("file_stem"),
            created_at=data.get("created_at", now_iso()),
        )

    def to_dict(self) -> Dict[str, Any]:
        """Serializes the figure for dcc.Store or JSON export."""
        return {
            "id": self.id,
            "dataset_key": self.dataset_key,
            "view_id": self.view_id,
            "filter_state": self.filter_state.to_dict(),
            "view_params": self.view_params,
            "label": self.label,
            "file_stem": self.file_stem,
            "created_at": self.created_at,
        }


# -------------------------------------------------------------------------
# Session Metadata (Overarching Container)
# -------------------------------------------------------------------------

@dataclass
class SessionMetadata:
    """
    Overarching container for a session, matching your metadata 6.json structure.
    """
    session_id: str
    schema_version: int
    app_version: str
    datasets_config_hash: str
    created_at: str = field(default_factory=now_iso)
    updated_at: str = field(default_factory=now_iso)
    figures: List[FigureMetadata] = field(default_factory=list)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> SessionMetadata:
        """Reconstructs a full session from a dictionary."""
        figures = [FigureMetadata.from_dict(f) for f in data.get("figures", [])]

        return cls(
            session_id=data.get("session_id", generate_session_id()),
            schema_version=data.get("schema_version", 1),
            app_version=data.get("app_version", "0.0.0-dev"),
            datasets_config_hash=data.get("datasets_config_hash", "unknown"),
            created_at=data.get("created_at", now_iso()),
            updated_at=data.get("updated_at", now_iso()),
            figures=figures,
        )

    def to_dict(self) -> Dict[str, Any]:
        """Serializes the full session to a dictionary."""
        return {
            "session_id": self.session_id,
            "schema_version": self.schema_version,
            "app_version": self.app_version,
            "datasets_config_hash": self.datasets_config_hash,
            "created_at": self.created_at,
            "updated_at": self.updated_at,
            "figures": [f.to_dict() for f in self.figures],
        }


# -------------------------------------------------------------------------
# Creation Helper
# -------------------------------------------------------------------------

def new_session_metadata(
        session_id: Optional[str] = None,
        app_version: str = "0.0.0-dev",
        datasets_config_hash: str = "unknown",
) -> SessionMetadata:
    """Create a fresh SessionMetadata container."""
    now = now_iso()
    return SessionMetadata(
        session_id=session_id or generate_session_id(),
        schema_version=1,
        app_version=app_version,
        datasets_config_hash=datasets_config_hash,
        created_at=now,
        updated_at=now,
        figures=[],
    )