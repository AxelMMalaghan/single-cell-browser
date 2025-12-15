from __future__ import annotations

from pathlib import Path
from typing import Dict

import logging
import plotly.graph_objs as go

from sc_browser.core import Dataset
from sc_browser.core.filter_state import FilterState
from sc_browser.core.view_registry import ViewRegistry
from sc_browser.metadata_io.model import FigureMetadata

logger = logging.getLogger(__name__)


class ExportService:
    """
    Responsible for:
    - turning figure metadata into a Plotly figure
    - writing figures out as image files
    """

    def __init__(
            self,
            *,
            datasets_by_key: Dict[str, Dataset],
            view_registry: ViewRegistry,
            output_root: Path,
    ) -> None:
        self._datasets_by_key = datasets_by_key
        self._view_registry = view_registry
        self._output_root = output_root


    def _get_dataset(self, key: str) -> Dataset:
        try:
            return self._datasets_by_key[key]
        except KeyError:
            raise KeyError(f"Dataset with key {key} not found")


    def _create_view(self, view_id: str, dataset: Dataset):
            return self._view_registry.create(view_id, dataset)

    def render_figure(self, metadata: FigureMetadata) -> go.Figure:
        logger.debug(
            "Rendering figure: figure_id=%s dataset_key=%s view_id=%s",
            getattr(metadata, "id", None),
            metadata.dataset_key,
            metadata.view_id,
        )

        ds = self._get_dataset(metadata.dataset_key)
        view = self._create_view(metadata.view_id, ds)

        filter_state = metadata.filter_state
        if isinstance(filter_state, dict):
            logger.debug("filter_state provided as dict; converting via FilterState.from_dict")
            filter_state = FilterState.from_dict(filter_state)

        try:
            logger.debug(
                "Filter summary: clusters=%d conditions=%d samples=%d cell_types=%d genes=%d",
                len(getattr(filter_state, "clusters", []) or []),
                len(getattr(filter_state, "conditions", []) or []),
                len(getattr(filter_state, "samples", []) or []),
                len(getattr(filter_state, "cell_types", []) or []),
                len(getattr(filter_state, "genes", []) or []),
            )
        except Exception:
            # Don’t let logging break rendering
            pass

        data = view.compute_data(filter_state)
        try:
            logger.debug("Computed data: rows=%s cols=%s", getattr(data, "shape", ("?", "?"))[0],
                         getattr(data, "shape", ("?", "?"))[1])
        except Exception:
            logger.debug("Computed data (shape unknown)")

        figure = view.render_figure(data, filter_state)
        logger.debug("Render complete: figure_id=%s", getattr(metadata, "id", None))
        return figure

    def export_single(self, metadata: FigureMetadata, session_id: str, image_format: str = "png") -> Path:
        file_stem = metadata.file_stem or f"{metadata.dataset_key}.{metadata.view_id}_{metadata.id}"

        session_dir = self._output_root / session_id
        out_path = session_dir / f"{file_stem}.{image_format}"

        logger.info(
            "Exporting figure: session_id=%s figure_id=%s dataset_key=%s view_id=%s format=%s out_path=%s",
            session_id,
            getattr(metadata, "id", None),
            metadata.dataset_key,
            metadata.view_id,
            image_format,
            out_path,
        )

        try:
            figure = self.render_figure(metadata)
            session_dir.mkdir(parents=True, exist_ok=True)
            figure.write_image(str(out_path))
        except Exception:
            logger.exception(
                "Export failed: session_id=%s figure_id=%s out_path=%s",
                session_id,
                getattr(metadata, "id", None),
                out_path,
            )
            raise

        logger.info("Export complete: out_path=%s", out_path)
        return out_path

