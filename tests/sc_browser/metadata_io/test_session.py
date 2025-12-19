from __future__ import annotations

from sc_browser.core.filter_state import FilterState
from sc_browser.core.metadata_model import (
    new_session_metadata,
    session_to_dict,
    session_from_dict,
    generate_session_id,
)
from sc_browser.core.metadata_model import FigureMetadata


def test_session_roundtrip_preserves_filter_state_type():
    session = new_session_metadata(
        session_id=generate_session_id(),
        app_version="0.0-test",
        datasets_config_hash="abc123",
    )

    state = FilterState(dataset_name="ds1", view_id="cluster")
    fig = FigureMetadata.from_runtime(
        figure_id="fig-0001",
        dataset_key="ds1",
        view_id="cluster",
        state=state,
        view_params={},
        label="Test",
        file_stem="fig-0001",
    )
    session.figures.append(fig)

    raw = session_to_dict(session)
    rebuilt = session_from_dict(raw)

    assert rebuilt is not None
    assert len(rebuilt.figures) == 1
    assert isinstance(rebuilt.figures[0].filter_state, FilterState), (
        "filter_state must be a FilterState, not a dict. "
        "Fix session_from_dict to call FilterState.from_dict(...)"
    )
