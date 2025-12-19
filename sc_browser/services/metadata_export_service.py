from __future__ import annotations

import io
import json
import logging
import zipfile
from typing import Dict, Any

import plotly.graph_objs as go

from sc_browser.core import Dataset
from sc_browser.core.filter_state import FilterState
from sc_browser.core.view_registry import ViewRegistry
from sc_browser.metadata_io.metadata_model import FigureMetadata, SessionMetadata

logger = logging.getLogger(__name__)


class ExportService:
    """
    Handles rendering figures and packaging them into a ZIP bundle.
    Stateless: renders images from metadata on-the-fly for export.
    """

    def __init__(
            self,
            *,
            datasets_by_key: Dict[str, Dataset],
            view_registry: ViewRegistry,
    ) -> None:
        self._datasets_by_key = datasets_by_key
        self._view_registry = view_registry

    def _get_dataset(self, key: str) -> Dataset:
        try:
            return self._datasets_by_key[key]
        except KeyError:
            # Fallback to name-based lookup if key is missing
            raise KeyError(f"Dataset with key {key} not found")

    def render_figure(self, metadata: FigureMetadata) -> go.Figure:
        """Renders a Plotly figure object from FigureMetadata."""
        ds = self._get_dataset(metadata.dataset_key)
        view = self._view_registry.create(metadata.view_id, ds)

        # Ensure filter_state is a FilterState object (fixes reconstruction bug)
        filter_state = metadata.filter_state
        if isinstance(filter_state, dict):
            filter_state = FilterState.from_dict(filter_state)

        data = view.compute_data(filter_state)
        return view.render_figure(data, filter_state)

    def create_session_zip(self, session_data: Dict[str, Any]) -> bytes:
        """
        Generates a ZIP bundle containing metadata.json and rendered PNGs.
        Stateless: iterates through the provided dictionary and renders on-the-fly.
        """
        buf = io.BytesIO()

        # Reconstruct SessionMetadata to validate and use helper methods
        session = SessionMetadata.from_dict(session_data)

        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
            # 1. Metadata - Use the new class method
            # FIX: Replaced session_to_dict(session) with session.to_dict()
            zf.writestr("metadata.json", json.dumps(session.to_dict(), indent=2).encode("utf-8"))

            # 2. Images - Render each figure to PNG and add to the ZIP
            added_count = 0
            for fig_meta in session.figures:
                try:
                    # Render the figure object
                    figure = self.render_figure(fig_meta)

                    # Convert to PNG bytes (requires kaleido or orca)
                    img_bytes = figure.to_image(format="png")

                    # Name the file using stem or fallback pattern
                    filename = fig_meta.file_stem or f"{fig_meta.dataset_key}.{fig_meta.view_id}_{fig_meta.id}"
                    zf.writestr(f"figures/{filename}.png", img_bytes)
                    added_count += 1
                except Exception as e:
                    logger.error(f"Failed to render figure {fig_meta.id} for export: {e}")

            if added_count == 0 and session.figures:
                zf.writestr("WARNING.txt", b"No images could be rendered for this session.")

        return buf.getvalue()