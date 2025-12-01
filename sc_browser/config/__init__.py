"""
Config package for sc_browser.

Responsible for:
- config models (GlobalConfig, DatasetConfig, etc.)
- config I/O helpers (load_datasets / load_global_config)
"""

from .model import GlobalConfig, DatasetConfig  # optional re-exports
from .io import load_global_config, load_datasets  # optional