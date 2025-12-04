from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any
from plotly.graph_objs import Figure

from .dataset import Dataset
from .state import FilterState, FilterProfile


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
    def render_figure(self, data: Any, state: FilterState) -> Figure:
        """
        Render the figure given the computed data
        :param data: the data provided by {@link compute_data()}
        :param state: the current state of the {@link FilterState} object - what filters the user has toggled for
        :return: the Plotly figure for these parameters
        """
        raise NotImplementedError()

