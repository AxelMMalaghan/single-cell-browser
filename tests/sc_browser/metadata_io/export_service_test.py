from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
import plotly.graph_objs as go

from sc_browser.core.filter_state import FilterState
from sc_browser.metadata_io.metadata_model import FigureMetadata
from sc_browser.services.metadata_export_service import ExportService


class _DummyView:
    def __init__(self, dataset):
        self.dataset = dataset
        self.compute_data_called = 0
        self.render_figure_called = 0
        self.last_filter_state = None
        self.last_data = None
        self.last_render_filter_state = None

    def compute_data(self, filter_state):
        self.compute_data_called += 1
        self.last_filter_state = filter_state
        return pd.DataFrame({"x": [1], "y": [2]})

    def render_figure(self, data, filter_state):
        self.render_figure_called += 1
        self.last_data = data
        self.last_render_filter_state = filter_state
        return go.Figure()

class _DummyRegistry:
    def __init__(self):
        self.last_view = None

    def create(self, view_id, dataset):
        self.last_view = _DummyView(dataset)
        return self.last_view


class _DummyDataset:
    def __init__(self, name: str = "ds1") -> None:
        self.name = name


@pytest.fixture()
def dataset_map():
    # ExportService uses dataset_key lookup; keep it simple
    return {"ds1": _DummyDataset("ds1")}


@pytest.fixture()
def state() -> FilterState:
    return FilterState(dataset_name="ds1", view_id="cluster")


def test_render_figure_calls_view_and_returns_plotly_figure(tmp_path: Path, dataset_map, state: FilterState):
    registry = _DummyRegistry()

    svc = ExportService(
        datasets_by_key=dataset_map,
        view_registry=registry,
        output_root=tmp_path,
    )

    md = FigureMetadata.from_runtime(
        figure_id="fig-0001",
        dataset_key="ds1",
        view_id="cluster",
        state=state,
        view_params={},
        label="Test",
        file_stem="fig-0001",
    )

    fig = svc.render_figure(md)

    expected_state = md.filter_state
    if isinstance(expected_state, dict):
        expected_state = FilterState.from_dict(expected_state)

    assert registry.last_view.compute_data_called == 1
    assert registry.last_view.render_figure_called == 1
    assert registry.last_view.last_filter_state == expected_state
    assert registry.last_view.last_render_filter_state == expected_state
    assert isinstance(registry.last_view.last_data, pd.DataFrame)

def test_write_image_is_mockable(tmp_path: Path, monkeypatch, dataset_map, state: FilterState):
    registry = _DummyRegistry()
    svc = ExportService(
        datasets_by_key=dataset_map,
        view_registry=registry,
        output_root=tmp_path,
    )

    md = FigureMetadata.from_runtime(
        figure_id="fig-0001",
        dataset_key="ds1",
        view_id="cluster",
        state=state,
        view_params={},
        label="Test",
        file_stem="fig-0001",
    )

    fig = svc.render_figure(md)

    called = {"path": None}

    def _fake_write_image(self, file, *args, **kwargs):
        called["path"] = Path(file)

    monkeypatch.setattr(go.Figure, "write_image", _fake_write_image, raising=True)

    out = tmp_path / "fig.png"
    fig.write_image(out)

    assert called["path"] == out