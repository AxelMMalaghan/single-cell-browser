from __future__ import annotations

import json
import logging
from typing import Any, Optional, Tuple

from sc_browser.metadata_io.io import normalise_session_dict
from sc_browser.metadata_io.metadata_model import (
    FigureMetadata,
    SessionMetadata,
    generate_figure_id,
    new_session_metadata,
    now_iso,
    session_to_dict,
)
from sc_browser.services.storage import StorageBackend

logger = logging.getLogger(__name__)


class SessionService:
    """
    Manages session lifecycle and persistence.
    Ensures that session metadata is backed up to storage on every modification.
    """

    def __init__(self, storage: StorageBackend):
        self.storage = storage

    def _persist(self, session: SessionMetadata) -> None:
        """Internal helper to save session to disk."""
        try:
            # We save under the session_id folder
            path = f"{session.session_id}/metadata.json"
            data = session_to_dict(session)
            json_bytes = json.dumps(data, indent=2).encode("utf-8")
            self.storage.write_bytes(path, json_bytes)
        except Exception:
            logger.exception("Failed to persist session %s", session.session_id)

    def persist_session(self, session: SessionMetadata) -> None:
        """Public method to force a save (e.g. after import/merge)."""
        self._persist(session)

    def load_session(self, session_id: str) -> Optional[SessionMetadata]:
        """Load session from storage if it exists."""
        path = f"{session_id}/metadata.json"
        if not self.storage.exists(path):
            return None
        try:
            data = json.loads(self.storage.read_bytes(path))
            return normalise_session_dict(data)
        except Exception:
            logger.exception("Failed to load session %s", session_id)
            return None

    def ensure_session(
            self,
            session: Optional[SessionMetadata],
            *,
            session_id: str,
            app_version: str = "0.0.0-dev",
            datasets_config_hash: str = "unknown",
    ) -> SessionMetadata:
        """
        Ensure a SessionMetadata exists. If creating new, persists it immediately.
        """
        if session is not None:
            return session

        # Try loading from disk first (recovery)
        loaded = self.load_session(session_id)
        if loaded:
            return loaded

        # Create fresh
        new_sess = new_session_metadata(
            session_id=session_id,
            app_version=app_version,
            datasets_config_hash=datasets_config_hash,
        )
        self._persist(new_sess)
        return new_sess


    def save_figure(
            self,
            session: SessionMetadata,
            *,
            active_figure_id: Optional[str],
            dataset_key: str,
            view_id: str,
            filter_state: dict[str, Any],
            label: Any,
    ) -> Tuple[SessionMetadata, str, bool]:

        label_clean = str(label).strip() or None if label else None

        # Intent: If the UI provides an ID, we try to update that specific figure.
        existing_idx = None
        if active_figure_id:
            existing_idx = next((i for i, f in enumerate(session.figures) if f.id == active_figure_id), None)

        is_overwrite = existing_idx is not None

        meta = FigureMetadata.from_runtime(
            figure_id=active_figure_id if is_overwrite else generate_figure_id(),
            dataset_key=dataset_key,
            view_id=view_id,
            state=filter_state,
            view_params={},
            label=label_clean,
        )

        if is_overwrite:
            session.figures[existing_idx] = meta
        else:
            session.figures.append(meta)

        session.updated_at = now_iso()
        self._persist(session)
        return session, meta.id, is_overwrite

    def delete_figure(self, session: SessionMetadata, *, figure_id: str) -> SessionMetadata:
        """
        Remove a figure from the session and persist to disk.
        """
        initial_len = len(session.figures)
        session.figures = [f for f in session.figures if f.id != figure_id]

        if len(session.figures) != initial_len:
            session.updated_at = now_iso()
            self._persist(session)

        return session

