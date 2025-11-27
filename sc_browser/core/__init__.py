"""
Core domain layer: dataset abstraction, filter state, view base class,
and the view registry
"""

from .dataset import Dataset
from .state import FilterState
from .view_base import BaseView
from .view_registry import ViewRegistry

__all__ = ["Dataset", "FilterState", "BaseView", "ViewRegistry"]