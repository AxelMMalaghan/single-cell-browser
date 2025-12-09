from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict

from .model import FigureMetadata, SessionMetadata, session_from_dict, new_session_metadata, _now_iso


def load_session_metadata_from_file(path: Path) -> SessionMetadata:
    """
    Loads the session metadata from a file path
    :param path: the file path
    :return: the session metadata
    """

    with path.open() as f:
        raw: Dict[str, Any] = json.load(f)

    return normalise_session_dict(raw)

def normalise_session_dict(raw: Dict[str, Any]) -> SessionMetadata:
    """
    Normalises the session metadata from a raw file
    :param raw: the raw metadata file
    :return: session metadata
    """
    # Case 1 - Already a full session
    if "session_id" in raw and "figures" in raw:
        session = session_from_dict(raw)

        if session is None:
            raise ValueError("Invalid session metadata JSON")
        return session


    # Case 2 - bare bundle:
    if "figures" in raw:
        figures = [
            FigureMetadata(
                id = f["id"],
                dataset_key=f["dataset_key"],
                view_id=f["view_id"],
                filter_state=f["filter_state"],
                view_params=f["view_params"],
                label=f["label"],
                file_stem=f["file_stem"]
            ) for f in raw["figures"]
        ]
    else:
        # Case 3 single-figure json
        figures = [
            FigureMetadata(
                id=["id"],
                dataset_key=["dataset_key"],
                view_id=["view_id"],
                filter_state=["filter_state"],
                view_params=["view_params"],
                label=["label"],
                file_stem=["file_stem"]
            )
        ]

    session = new_session_metadata(
        session_id=f"import-{_now_iso()}",
        app_version="import-unknown",
        datasets_config_hash="import-unknown",
    )
    session.figures.extend(figures)
    return session
