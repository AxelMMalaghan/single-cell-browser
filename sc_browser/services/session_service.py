
from __future__ import annotations

import json
from dataclasses import replace
from typing import Any, Iterable, Optional, Tuple

from sc_browser.metadata_io.model import (
    FigureMetadata,
    SessionMetadata,
    generate_figure_id,
    new_session_metadata,
    now_iso,
)


def _clean_label(label: Any) -> Optional[str]:
    if label is None:
        return None
    s = str(label).strip()
    return s or None


def _figure_content_key(fig: FigureMetadata) -> tuple:
    """
    Content-based key for dedupe (NOT using fig.id).
    """
    return (
        getattr(fig, "dataset_key", None),
        getattr(fig, "view_id", None),
        json.dumps(getattr(fig, "filter_state", None) or {}, sort_keys=True),
        getattr(fig, "label", None),
    )


def ensure_session(
    session: Optional[SessionMetadata],
    *,
    session_id: str,
    app_version: str = "0.0.0-dev",
    datasets_config_hash: str = "unknown",
) -> SessionMetadata:
    """
    Ensure a SessionMetadata exists.
    """
    if session is not None:
        return session
    return new_session_metadata(
        session_id=session_id,
        app_version=app_version,
        datasets_config_hash=datasets_config_hash,
    )


def save_figure(
    session: SessionMetadata,
    *,
    active_figure_id: Optional[str],
    dataset_key: str,
    view_id: str,
    filter_state: dict[str, Any],
    label: Any,
) -> Tuple[SessionMetadata, str, bool]:
    """
    Save current view state into session.figures.

    Overwrite rule matches your current behavior:
    - If active_figure_id exists AND label matches existing label => overwrite that id
    - Else create a new id

    Returns: (updated_session, new_active_figure_id, is_overwrite)
    """
    label_clean = _clean_label(label)

    existing_fig: Optional[FigureMetadata] = None
    if active_figure_id:
        existing_fig = next((f for f in session.figures if f.id == active_figure_id), None)

    same_label_as_existing = (
        existing_fig is not None and (label_clean or "") == (existing_fig.label or "")
    )

    if existing_fig is not None and same_label_as_existing:
        figure_id = existing_fig.id
        is_overwrite = True
    else:
        figure_id = generate_figure_id(session)
        is_overwrite = False

    meta = FigureMetadata.from_runtime(
        figure_id=figure_id,
        dataset_key=dataset_key,
        view_id=view_id,
        state=filter_state,
        view_params={},
        label=label_clean,
        file_stem=None,
    )

    if is_overwrite and existing_fig is not None:
        # replace in-place by id
        for idx, f in enumerate(session.figures):
            if f.id == existing_fig.id:
                session.figures[idx] = meta
                break
    else:
        session.figures.append(meta)

    session.updated_at = now_iso()
    return session, figure_id, is_overwrite


def delete_figure(session: SessionMetadata, *, figure_id: str) -> SessionMetadata:
    """
    Remove a figure from the session by id.
    """
    session.figures = [f for f in session.figures if f.id != figure_id]
    session.updated_at = now_iso()
    return session


def merge_sessions(
    existing: Optional[SessionMetadata],
    imported: SessionMetadata,
    *,
    max_figures: int = 100,
) -> Tuple[SessionMetadata, int]:
    """
    Merge existing + imported figures, dedupe by content, then re-id sequentially.

    Returns: (merged_session, n_imported_effective_estimate)
    """
    # Defensive truncate on input
    imported_figs = list(imported.figures[:max_figures])

    if existing is None:
        merged = imported
        existing_figs: list[FigureMetadata] = []
    else:
        merged = existing
        existing_figs = list(existing.figures)

    all_figs = existing_figs + imported_figs

    seen: set[tuple] = set()
    unique: list[FigureMetadata] = []
    for fig in all_figs:
        k = _figure_content_key(fig)
        if k in seen:
            continue
        seen.add(k)
        unique.append(fig)

    # Re-id sequentially without mutating frozen dataclasses
    normalised: list[FigureMetadata] = []
    for idx, fig in enumerate(unique[:max_figures], start=1):
        normalised.append(replace(fig, id=f"fig-{idx:04d}"))

    merged.figures = normalised
    merged.updated_at = now_iso()

    # best-effort: how many imported survived
    # (not perfect, but good enough for UX banner)
    n_imported_effective = min(len(imported_figs), max_figures)

    return merged, n_imported_effective
