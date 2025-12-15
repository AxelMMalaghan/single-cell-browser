from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List

from .model import FigureMetadata, SessionMetadata, session_from_dict, new_session_metadata, now_iso


def load_session_metadata_from_file(path: Path) -> SessionMetadata:
    """
    Loads the session metadata from a file path.
    """
    with path.open() as f:
        raw: Dict[str, Any] = json.load(f)

    return normalise_session_dict(raw)


def _renumber_figures_sequentially(figures: List[FigureMetadata]) -> None:
    """
    Enforce unique, purely sequential ids: fig-0001, fig-0002, ...

    This intentionally discards imported ids to avoid collisions.
    """
    for i, fig in enumerate(figures, start=1):
        fig.id = f"fig-{i:04d}"


def normalise_session_dict(raw: Dict[str, Any]) -> SessionMetadata:
    """
    Normalises the session metadata from a raw file.
    Supports:
      1) full session JSON: {"session_id": ..., "figures": [...]}
      2) bare bundle: {"figures": [...]}
      3) single-figure JSON: {"id": ..., "dataset_key": ..., ...}
    """
    # Case 1 - Already a full session
    if "session_id" in raw and "figures" in raw:
        session = session_from_dict(raw)
        if session is None:
            raise ValueError("Invalid session metadata JSON")

        _renumber_figures_sequentially(session.figures)
        return session

    figures: List[FigureMetadata] = []

    # Case 2 - bare bundle:
    if "figures" in raw and isinstance(raw["figures"], list):
        for f in raw["figures"]:
            if not isinstance(f, dict):
                continue
            figures.append(
                FigureMetadata(
                    id=str(f.get("id", "")).strip() or f"fig-{len(figures)+1:04d}",  # overwritten below
                    dataset_key=str(f.get("dataset_key", "")).strip(),
                    view_id=str(f.get("view_id", "")).strip(),
                    filter_state=f.get("filter_state") or {},
                    view_params=f.get("view_params") or {},
                    label=f.get("label"),
                    file_stem=f.get("file_stem"),
                    created_at=f.get("created_at", now_iso()),
                )
            )
    else:
        # Case 3 - single-figure json
        figures.append(
            FigureMetadata(
                id=str(raw.get("id", "")).strip() or "fig-0001",  # overwritten below
                dataset_key=str(raw.get("dataset_key", "")).strip(),
                view_id=str(raw.get("view_id", "")).strip(),
                filter_state=raw.get("filter_state") or {},
                view_params=raw.get("view_params") or {},
                label=raw.get("label"),
                file_stem=raw.get("file_stem"),
                created_at=raw.get("created_at", now_iso()),
            )
        )

    _renumber_figures_sequentially(figures)

    session = new_session_metadata(
        session_id=f"import-{now_iso()}",
        app_version="import-unknown",
        datasets_config_hash="import-unknown",
    )
    session.figures.extend(figures)
    return session
