from __future__ import annotations

import logging
import os
from typing import Optional

from pythonjsonlogger import jsonlogger

def configure_logging(
        level: int = logging.INFO,
        force_format: Optional[str] = None,
) -> None:
    """
    Configure root logger for the app

    Modes:
    - JSON (default) in prod
    - plain text (dev mode)

    Selection Order:
        1) force_format argument ("json or "plain") if provided
        2) env var SC_BROWSER_LOG_FORMAT
        3) default = "json"
    """

    if force_format is not None:
        format_mode = force_format
    else:
        format_mode = os.getenv("SC_BROWSER_LOG_FORMAT", "json").lower()

    logger = logging.getLogger()
    logger.setLevel(level)

    handler = logging.StreamHandler()

    if format_mode == "plain":
        # Test logs (dev mode)
        formatter = logging.Formatter(
            "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
        )
    else:
        formatter = jsonlogger.JsonFormatter(
            "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
        )

    handler.setFormatter(formatter)

    # Replace any existing handlers to avoid duplicate logs
    logger.handlers.clear()
    logger.addHandler(handler)

