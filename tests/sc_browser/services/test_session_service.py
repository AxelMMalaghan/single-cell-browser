from __future__ import annotations

from sc_browser.services.session_service import (
    ensure_session,
    save_figure,
    delete_figure,
    merge_sessions,
)
from sc_browser.core.metadata_model import new_session_metadata


def _new_session(session_id: str = "s1"):
    return new_session_metadata(
        session_id=session_id,
        app_version="test",
        datasets_config_hash="hash",
    )


def test_save_figure_creates_new_when_no_active():
    session = _new_session()
    session, fig_id, is_overwrite = save_figure(
        session,
        active_figure_id=None,
        dataset_key="ds1",
        view_id="heatmap",
        filter_state={"dataset_name": "X", "view_id": "heatmap"},
        label="My Figure",
    )
    assert is_overwrite is False
    assert fig_id.startswith("fig-")
    assert len(session.figures) == 1
    assert session.figures[0].id == fig_id
    assert session.figures[0].label == "My Figure"


def test_save_figure_overwrites_when_active_and_same_label():
    session = _new_session()

    session, fig_id1, _ = save_figure(
        session,
        active_figure_id=None,
        dataset_key="ds1",
        view_id="umap",
        filter_state={"a": 1},
        label="A",
    )

    # same active + same label => overwrite same id (no new figure)
    session, fig_id2, is_overwrite = save_figure(
        session,
        active_figure_id=fig_id1,
        dataset_key="ds1",
        view_id="umap",
        filter_state={"a": 2},
        label="A",
    )

    assert is_overwrite is True
    assert fig_id2 == fig_id1
    assert len(session.figures) == 1
    assert session.figures[0].id == fig_id1
    assert session.figures[0].filter_state == {"a": 2}


def test_delete_figure_removes_by_id():
    session = _new_session()
    session, fig_id, _ = save_figure(
        session,
        active_figure_id=None,
        dataset_key="ds1",
        view_id="umap",
        filter_state={"a": 1},
        label="A",
    )
    assert len(session.figures) == 1

    session = delete_figure(session, figure_id=fig_id)
    assert len(session.figures) == 0


def test_merge_sessions_dedupes_by_content_and_reids():
    existing = _new_session("s-existing")
    imported = _new_session("s-imported")

    existing, _, _ = save_figure(
        existing,
        active_figure_id=None,
        dataset_key="ds1",
        view_id="umap",
        filter_state={"x": 1},
        label="A",
    )

    # duplicate content (same dataset/view/state/label)
    imported, _, _ = save_figure(
        imported,
        active_figure_id=None,
        dataset_key="ds1",
        view_id="umap",
        filter_state={"x": 1},
        label="A",
    )

    # unique content
    imported, _, _ = save_figure(
        imported,
        active_figure_id=None,
        dataset_key="ds1",
        view_id="umap",
        filter_state={"x": 2},
        label="B",
    )

    merged, n_imported_eff = merge_sessions(existing, imported, max_figures=100)
    assert n_imported_eff == 2
    assert len(merged.figures) == 2

    # ids are sequential and normalised
    assert merged.figures[0].id == "fig-0001"
    assert merged.figures[1].id == "fig-0002"
