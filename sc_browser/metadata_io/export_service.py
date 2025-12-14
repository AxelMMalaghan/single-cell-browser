from __future__ import annotations

from pathlib import Path
from typing import Dict

import plotly.graph_objs as go

from sc_browser.core import Dataset
from sc_browser.core.filter_state import FilterState
from sc_browser.core.view_registry import ViewRegistry
from sc_browser.metadata_io.model import FigureMetadata


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
        """
        Render figure based on figure metadata
        :param metadata: the figure metadata instance
        :return: a plotly figure
        """
        ds = self._get_dataset(metadata.dataset_key)
        view = self._create_view(metadata.view_id, ds)

        filter_state = metadata.filter_state
        if isinstance(filter_state, dict):
            filter_state = FilterState.from_dict(filter_state)

        data = view.compute_data(filter_state)
        figure = view.render_figure(data, filter_state)

        return figure


    def export_single(
            self,
            metadata: FigureMetadata,
            session_id: str,
            image_format: str = "png",
    ) -> Path:
        """
        Render figure based on figure metadata
        :param metadata: the figure metadata instance
        :param session_id: the session id
        :param image_format: the image format
        :return: the path to the exported image
        """

        figure = self.render_figure(metadata)

        session_dir = self._output_root / session_id
        session_dir.mkdir(parents=True, exist_ok=True)

        file_stem = metadata.file_stem or f"{metadata.dataset_key}.{metadata.view_id}_{metadata.id}"
        out_path = session_dir / f"{file_stem}.{image_format}"

        figure.write_image(str(out_path))
        return out_path
