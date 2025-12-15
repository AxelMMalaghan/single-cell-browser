from __future__ import annotations
import logging
from typing import Optional

from sc_browser.core.filter_state import FilterState

logger = logging.getLogger(__name__)

def try_parse_filter_state(data: object) -> Optional[FilterState]:
    if not isinstance(data, dict) or not data:
        return None
    try:
        return FilterState.from_dict(data)
    except Exception:
        logger.exception("Invalid filter-state: %r", data)
        return None


def safe_filter_state(data) -> Optional[FilterState]:
    if not isinstance(data, dict):
        return None
    try:
        return FilterState.from_dict(data)
    except Exception:
        return None