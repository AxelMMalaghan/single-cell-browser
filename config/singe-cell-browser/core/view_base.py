from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Any

import pandas as pd
from plotly.graph_objs import Figure  # type: ignore

from .dataset import Dataset
from .state import FilterState


class base_view(ABC):
    """
    Base class for all plots

    Each subclass must define:
    - id
    - label
    - compute_data()
    - render_figure()
    """

    id: str = None
    label: str = None

    def __init__(self, dataset: Dataset):
        self.dataset = dataset


    @abstractmethod
    def compute_data(self, state: FilterState) -> Any:
        raise NotImplementedError()

    @abstractmethod
    def render_figure(self, data: Any, state: FilterState) -> Figure:
        raise NotImplementedError()

