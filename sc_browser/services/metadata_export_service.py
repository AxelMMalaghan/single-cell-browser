from __future__ import annotations

import io
import json
import logging
import zipfile
from typing import Dict

import plotly.graph_objs as go

from sc_browser.core import Dataset
from sc_browser.core.filter_state import FilterState
from sc_browser.core.view_registry import ViewRegistry
from sc_browser.metadata_io.metadata_model import FigureMetadata, SessionMetadata, session_to_dict
from sc_browser.services.storage import StorageBackend

logger = logging.getLogger(__name__)


class ExportService:
    """
    Handles rendering and persistence of figures/sessions using an abstract storage backend.
    """

    def __init__(
            self,
            *,
            datasets_by_key: Dict[str, Dataset],
            view_registry: ViewRegistry,
            storage: StorageBackend,
    ) -> None:
        self._datasets_by_key = datasets_by_key
        self._view_registry = view_registry
        self._storage = storage

    def _get_dataset(self, key: str) -> Dataset:
        try:
            return self._datasets_by_key[key]
        except KeyError:
            raise KeyError(f"Dataset with key {key} not found")

    def render_figure(self, metadata: FigureMetadata) -> go.Figure:
        # ... (Existing render logic remains unchanged) ...
        ds = self._get_dataset(metadata.dataset_key)
        view = self._view_registry.create(metadata.view_id, ds)

        filter_state = metadata.filter_state
        if isinstance(filter_state, dict):
            filter_state = FilterState.from_dict(filter_state)

        data = view.compute_data(filter_state)
        return view.render_figure(data, filter_state)

    def export_single(self, metadata: FigureMetadata, session_id: str, image_format: str = "png") -> str:
        """
        Renders and saves a figure image to storage. Returns the relative path.
        """
        file_stem = metadata.file_stem or f"{metadata.dataset_key}.{metadata.view_id}_{metadata.id}"
        # We enforce a folder structure: session_id/filename
        rel_path = f"{session_id}/{file_stem}.{image_format}"

        logger.info("Exporting figure to storage: %s", rel_path)

        try:
            figure = self.render_figure(metadata)

            # Plotly write_image writes to bytes if file is str/Path,
            # but we need to bridge to our storage backend.
            # We use to_image() to get bytes directly.
            img_bytes = figure.to_image(format=image_format)

            self._storage.write_bytes(rel_path, img_bytes)

        except Exception:
            logger.exception("Export failed for %s", rel_path)
            raise

        return rel_path

    def create_session_zip(self, session: SessionMetadata) -> bytes:
        """
        Generates the export ZIP bundle for a session.
        Encapsulates file gathering logic so UI doesn't need to touch storage.
        """
        buf = io.BytesIO()
        session_id = session.session_id

        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
            # 1. Metadata
            meta_dict = session_to_dict(session)
            zf.writestr("metadata.json", json.dumps(meta_dict, indent=2).encode("utf-8"))

            # 2. Images
            # Logic: Scan storage for files belonging to this session
            # We look in the session_id directory

            # Get list of files in this session folder
            available_files = self._storage.list_files(prefix=session_id, suffix=".png")

            # Build set of valid stems from current session metadata (ZOMBIE FIX logic)
            valid_stems = set()
            for fig in session.figures:
                if fig.file_stem:
                    valid_stems.add(fig.file_stem)
                else:
                    valid_stems.add(f"{fig.dataset_key}.{fig.view_id}_{fig.id}")

            added_count = 0
            for file_path in available_files:
                # file_path is like "session-123/dataset.view_fig-001.png"
                # extract stem
                filename = file_path.split("/")[-1]
                stem = filename.rsplit(".", 1)[0]

                if stem in valid_stems:
                    img_data = self._storage.read_bytes(file_path)
                    zf.writestr(f"figures/{filename}", img_data)
                    added_count += 1

            if added_count == 0 and session.figures:
                zf.writestr("WARNING.txt", b"No matching images found in storage for this session.")

        return buf.getvalue()