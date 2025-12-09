from __future__ import annotations

from typing import Dict, Any, Tuple
from sc_browser.core.filter_state import FilterState
from sc_browser.metadata_io.model import FigureMetadata

def figure_to_runtime(metadata: FigureMetadata) -> Tuple[str, str, Dict[str, Any], Dict[str, Any]]:
    """
    Converts FigureMetadata back into runtime objects
    :param metadata: the figure metadata
    :return: tuple of (dataset_key, view_id, filter_state, view_params)
    """

    return (
        metadata.dataset_key,
        metadata.view_id,
        metadata.filter_state,
        metadata.view_params,
    )


def figure_to_filter_state(metadata: FigureMetadata) -> FilterState:
    """
    Converts FigureMetadata back into FilterState object to be parsed into a view
    :param metadata: the figure metadata
    :return: the filter state (user filter toggles)
    """
    return FilterState(**metadata.filter_state)