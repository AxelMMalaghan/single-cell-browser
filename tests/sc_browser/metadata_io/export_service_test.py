from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from sc_browser.metadata_io.export_service import ExportService
from sc_browser.core.filter_state import FilterState


def test_render_figure_builds_filter_state_from_dict_and_calls_view_pipeline():
    ds = object()

    view = MagicMock()
    view.compute_data.return_value = {"some": "data"}
    fig = MagicMock()
    view.render_figure.return_value = fig

    view_registry = MagicMock()
    view_registry.create.return_value = view

    svc = ExportService(
        datasets_by_key={"ds1": ds},
        view_registry=view_registry,
        output_root=None,  # not used here
    )

    fs = FilterState(
        dataset_name="ds1",
        view_id="v1",
        genes=["MS4A1"],
        clusters=["T"],
        conditions=["treated"],
        merge_genes=True,
        split_by_condition=False,
        is_3d=True,
        dim_x=1,
        dim_y=2,
        dim_z=3,
        color_scale="viridis",
    )

    metadata = SimpleNamespace(
        dataset_key="ds1",
        view_id="v1",
        filter_state=fs.to_dict(),
        file_stem=None,
        id="fig-0001",
    )

    out = svc.render_figure(metadata)

    view_registry.create.assert_called_once_with("v1", ds)

    # Assert FilterState instance & contents passed through
    view.compute_data.assert_called_once()
    (fs_arg,) = view.compute_data.call_args.args
    assert isinstance(fs_arg, FilterState)
    assert fs_arg.to_dict() == fs.to_dict()

    view.render_figure.assert_called_once_with({"some": "data"}, fs_arg)
    assert out is fig


def test_render_figure_missing_dataset_key_raises_keyerror():
    view_registry = MagicMock()
    svc = ExportService(datasets_by_key={}, view_registry=view_registry, output_root=None)

    fs = FilterState(dataset_name="nope", view_id="v1")
    metadata = SimpleNamespace(
        dataset_key="nope",
        view_id="v1",
        filter_state=fs.to_dict(),
        file_stem=None,
        id="fig-0001",
    )

    with pytest.raises(KeyError, match="Dataset with key nope not found"):
        svc.render_figure(metadata)

    view_registry.create.assert_not_called()


def test_export_single_writes_image_default_stem(tmp_path):
    ds = object()

    fig = MagicMock()
    view = MagicMock()
    view.compute_data.return_value = {"d": 1}
    view.render_figure.return_value = fig

    view_registry = MagicMock()
    view_registry.create.return_value = view

    svc = ExportService(
        datasets_by_key={"ds1": ds},
        view_registry=view_registry,
        output_root=tmp_path,
    )

    fs = FilterState(dataset_name="ds1", view_id="umap")
    metadata = SimpleNamespace(
        dataset_key="ds1",
        view_id="umap",
        filter_state=fs.to_dict(),
        file_stem=None,   # triggers default naming
        id="fig-0007",
    )

    out_path = svc.export_single(metadata, session_id="session-abc", image_format="png")

    expected_dir = tmp_path / "session-abc"
    assert expected_dir.is_dir()

    expected_stem = "ds1.umap_fig-0007"
    assert out_path == expected_dir / f"{expected_stem}.png"

    fig.write_image.assert_called_once_with(str(out_path))


def test_export_single_uses_custom_stem_and_format(tmp_path):
    ds = object()

    fig = MagicMock()
    view = MagicMock()
    view.compute_data.return_value = {"d": 1}
    view.render_figure.return_value = fig

    view_registry = MagicMock()
    view_registry.create.return_value = view

    svc = ExportService(
        datasets_by_key={"ds1": ds},
        view_registry=view_registry,
        output_root=tmp_path,
    )

    fs = FilterState(dataset_name="ds1", view_id="heatmap")
    metadata = SimpleNamespace(
        dataset_key="ds1",
        view_id="heatmap",
        filter_state=fs.to_dict(),
        file_stem="my_custom_name",
        id="fig-0002",
    )

    out_path = svc.export_single(metadata, session_id="session-xyz", image_format="jpg")

    expected = tmp_path / "session-xyz" / "my_custom_name.jpg"
    assert out_path == expected

    fig.write_image.assert_called_once_with(str(expected))