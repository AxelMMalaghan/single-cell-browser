from __future__ import annotations

import uuid
from dataclasses import dataclass, field
from typing import Any, Optional, Dict, List
from datetime import datetime, timezone

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

def generate_figure_id() -> str:
    """
    Generate a globally unique, immutable ID for figures
    """
    return f"fig-{uuid.uuid4().hex[:12]}"

def generate_session_id() -> str:
    """
    Generate a globally unique, immutable ID for sessions
    """
    return f"fig-{uuid.uuid4().hex[:8]}"

def now_iso() -> str:
    """
    Return a current UTC timestamp in ISO-8601 format.
    """
    return datetime.now(timezone.utc).isoformat()

# -------------------------------------------------------------------------
# Per-figure metadata
# -------------------------------------------------------------------------

@dataclass
class FigureMetadata:
    """
    Immutable description of a single saved figure.

    - id: stable identifier within a session (e.g. "fig-0001")
    - dataset_key: key/name used to look up Dataset from AppContext
    - view_id: BaseView.id ("cluster", "expression", etc.)
    - filter_state: serialized FilterState (as dict)
    - view_params: extra per-view knobs (embedding_key, color_by, etc.)
    - label: optional human-readable label shown in the Reports UI
    - file_stem: base filename used when exporting image/JSON
    - created_at: ISO8601 timestamp (UTC) when this figure was first saved
    """

    id: str
    dataset_key: str
    view_id: str
    filter_state: Dict[str, Any]
    view_params: Dict[str, Any] = field(default_factory=dict)

    label: Optional[str] = None
    file_stem: Optional[str] = None
    created_at: str = field(default_factory=now_iso)

    @classmethod
    def from_runtime(
        cls,
        *,
        figure_id: str,
        dataset_key: str,
        view_id: str,
        state: Dict[str, Any],
        view_params: Dict[str, Any],
        label: Optional[str] = None,
        file_stem: Optional[str] = None,
    ) -> "FigureMetadata":
        """
        Build a FigureMetadata from the current runtime state.
        """
        return cls(
            id=figure_id,
            dataset_key=dataset_key,
            view_id=view_id,
            filter_state=state or {},
            view_params=view_params or {},
            label=label,
            file_stem=file_stem,
        )


# -------------------------------------------------------------------------
# Session-level metadata (bundle of figures)
# -------------------------------------------------------------------------

@dataclass
class SessionMetadata:
    """
    Describes a reporting/metadata_io session.

    - session_id: stable identifier (can be random UUID or timestamp-based)
    - schema_version: version of this metadata schema (start at 1)
    - app_version: app version string
    - datasets_config_hash: hash of the datasets config used for this session
    - created_at: when the session was first created
    - updated_at: last time the session was modified (figure added/removed, etc.)
    - figures: list of FigureMetadata entries
    """

    session_id: str
    schema_version: int
    app_version: str
    datasets_config_hash: str

    created_at: str = field(default_factory=now_iso)
    updated_at: str = field(default_factory=now_iso)

    figures: List[FigureMetadata] = field(default_factory=list)




# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

def new_session_metadata(
    *,
    session_id: str,
    app_version: str,
    datasets_config_hash: str,
    schema_version: int = 1,
) -> SessionMetadata:
    """
    Create a fresh SessionMetadata with no figures.
    """
    now = now_iso()
    return SessionMetadata(
        session_id=session_id,
        schema_version=schema_version,
        app_version=app_version,
        datasets_config_hash=datasets_config_hash,
        created_at=now,
        updated_at=now,
        figures=[],
    )


def touch_session(session: SessionMetadata) -> None:
    """
    Update the 'updated_at' timestamp after mutating the session.
    """
    session.updated_at = now_iso()


def session_to_dict(session: SessionMetadata) -> Dict[str, Any]:
    return {
        "session_id": session.session_id,
        "schema_version": session.schema_version,
        "app_version": session.app_version,
        "datasets_config_hash": session.datasets_config_hash,
        "created_at": session.created_at,
        "updated_at": session.updated_at,
        "figures": [
            {
                "id": f.id,
                "dataset_key": f.dataset_key,
                "view_id": f.view_id,
                "filter_state": f.filter_state or {},
                "view_params": f.view_params or {},
                "label": f.label,
                "file_stem": f.file_stem,
                "created_at": f.created_at,
            }
            for f in session.figures
        ],
    }


def session_from_dict(
    data: Dict[str, Any] | None,
) -> SessionMetadata | None:
    """
    Rebuild SessionMetadata from a dict produced by session_to_dict.
    Returns None if data is None.
    """
    if data is None:
        return None

    figures = [
        FigureMetadata(
            id=f["id"],
            dataset_key=f["dataset_key"],
            view_id=f["view_id"],
            filter_state=f.get("filter_state") or {},
            view_params=f.get("view_params", {}),
            label=f.get("label"),
            file_stem=f.get("file_stem"),
            created_at=f.get("created_at", now_iso()),
        )
        for f in data.get("figures", [])
    ]

    return SessionMetadata(
        session_id=data["session_id"],
        schema_version=data["schema_version"],
        app_version=data["app_version"],
        datasets_config_hash=data["datasets_config_hash"],
        created_at=data.get("created_at", now_iso()),
        updated_at=data.get("updated_at", now_iso()),
        figures=figures,
    )


