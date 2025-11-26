from __future__ import annotations

from abc import ABC
from typing import Dict, List, Type

from .dataset import Dataset
from .view_base import BaseView

class ViewRegistry:
    """
    Registry for view classes so the app can build tabs dynamically
    """
    def __init__(self):
        self._views: Dict[str, Type[BaseView]] = {}


    def register(self, view_cls: Type[BaseView]) -> None:

        if not issubclass(view_cls, BaseView):
            raise TypeError("view_cls must be a subclass of BaseView")

        if view_cls.id in self._views:
            raise ValueError(f"View '{view_cls.id}' already registered")

        self._views[view_cls.id] = view_cls


    def create(self, view_id: str, dataset: Dataset) -> BaseView:

        try:
            cls = self._views[view_id]
        except KeyError:
            raise KeyError(f"View '{view_id}' not found")
        return cls(dataset)

    def all_classes(self) -> List[Type[BaseView]]:
        return list(self._views.values())
