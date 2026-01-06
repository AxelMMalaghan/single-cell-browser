from __future__ import annotations
from typing import Dict, List, Type

from .dataset import Dataset
from .base_view import BaseView


class ViewRegistry:
    """
    Registry for view classes so the app can build tabs dynamically

    Purpose:
    - Decouples UI/Dash layer from hardcoded view implementations by exposing {@link create(view_id, dataset)}
    - Enables dynamic construction of navigation (tabs, dropdowns, etc.) based on the registered views rather than hardcoded lists

    Design Notes:
    - Stores the subclasses of {@link BaseView}, not instances, so that each view can be instantiated on demand
    - Enforces variants:
        * only {@link BaseView} subclasses can be registered
        * each view 'id' is unique across the registry
    - Centralising creation logic makes it easier to extend for other behaviour later on (logging, DI) without changing calls
    """

    def __init__(self):
        self._views: Dict[str, Type[BaseView]] = {}

    def register(self, view_cls: Type[BaseView]) -> None:
        """
        Register a {@link BaseView} with the registry

        The class must:
        - be a subclass of {@link BaseView}
        - define a unique 'id' attribute

        :param view_cls: the subclass of {@link BaseView}

        Raises:
            TypeError: if view_cls is not a subclass of {@link BaseView}
            ValueError: if a view with same 'id' already exists
        """

        if not issubclass(view_cls, BaseView):
            raise TypeError(f"View '{view_cls.id}' must be a subclass of BaseView")

        if view_cls.id in self._views:
            raise ValueError(f"View '{view_cls.id}' already registered")

        self._views[view_cls.id] = view_cls

    def create(self, view_id: str, dataset: Dataset) -> BaseView:
        """
        Instantiate a view for the given view_id. Callers only need to know the 'id' of the view and pass in a dataset
        instance
        :param view_id: the id of the view
        :param dataset: the given dataset
        :return cls(): the instantiated view

        Raises:
            KeyError: if no view with the given id exists in the registry
        """
        try:
            cls = self._views[view_id]
        except KeyError:
            raise KeyError(f"View '{view_id}' not found")
        return cls(dataset)

    def all_classes(self) -> List[Type[BaseView]]:
        """
        Used at UI layer to build navigation elements - tab bars or dropdowns. Keeps UI fully driven by the registry.
        :return list: a list of view id's in the registry
        """
        return list(self._views.values())
