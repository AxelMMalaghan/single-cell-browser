from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any

import plotly.graph_objs as go

from .dataset import Dataset
from .filter_state import FilterState, FilterProfile


class BaseView(ABC):
    """
    Abstract base class for all plot views.

    Defines the contract that every view in the app must follow
    - expose an 'id' - used internally
    - expose a 'label' - used for UI/human-readable applications
    - implement 'compute_data' - used to compute the data given the current FilterState
    - implement 'render_figure' - used to render the figure using Plotly
    """

    id: str = None
    label: str = None
    filter_profile = FilterProfile()

    def __init__(self, dataset: Dataset):
        self.dataset = dataset


    @abstractmethod
    def compute_data(self, state: FilterState) -> Any:
        """
        Compute the data given the current FilterState
        :param state: the current state of the {@link FilterState} object - what filters the user has toggled for
        :return: data: a dataframe containing the data as per the FilterState
        """
        raise NotImplementedError()

    @abstractmethod
    def render_figure(self, data: Any, state: FilterState) -> go.Figure:
        """
        Render the figure given the computed data
        :param data: the data provided by {@link compute_data()}
        :param state: the current state of the {@link FilterState} object - what filters the user has toggled for
        :return: the Plotly figure for these parameters
        """
        raise NotImplementedError()


    # ------------------------------------------------------------------
    # Common helpers for all views
    # ------------------------------------------------------------------
    def filtered_dataset(self, state: FilterState) -> Dataset:
        """
        Return this view's Dataset filtered according to the given FilterState.

        All views should call this instead of touching .subset(...) directly,
        so if we ever need to change the filtering behaviour, we do it in one place.
        """
        return self.dataset.subset_for_state(state)

    @staticmethod
    def empty_figure(message: str) -> go.Figure:
        """
        Standardised 'no data' figure used by all views.
        """
        fig = go.Figure()
        fig.update_layout(
            title=message,
            xaxis={"visible": False},
            yaxis={"visible": False},
        )
        return fig