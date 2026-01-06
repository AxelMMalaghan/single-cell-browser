from __future__ import annotations

import logging
import os
import sys

from pythonjsonlogger.jsonlogger import JsonFormatter


def _env_bool(name: str, default: str = "0") -> bool:
    return os.getenv(name, default).strip().lower() in {"1", "true", "yes", "y", "on"}


def configure_logging() -> None:
    """
    Configure app logging for both dev and prod.

    Env vars:
      - LOG_LEVEL: DEBUG/INFO/WARNING/ERROR (default INFO, DEBUG if DEBUG=1)
      - LOG_FORMAT: "json" or "text" (default json in prod, text in dev)
      - DEBUG: if 1/true enables debug defaults
    """

    debug = _env_bool("DEBUG", "0")
    level_name = os.getenv("LOG_LEVEL", "DEBUG" if debug else "INFO").upper()
    level = getattr(logging, level_name, logging.INFO)

    # Heuristic: Docker/k8s/prod -> json, local dev -> text
    default_format = "text" if debug else "json"
    fmt = os.getenv("LOG_FORMAT", default_format).lower()

    root = logging.getLogger()
    root.setLevel(level)

    # Clear existing handlers to avoid duplicate logs (esp. with Dash reload / tests)
    root.handlers.clear()

    handler = logging.StreamHandler(sys.stdout)

    if fmt == "json":
        formatter = JsonFormatter(
            "%(asctime)s %(levelname)s %(name)s %(message)s"
        )
    else:
        formatter = logging.Formatter(
            fmt="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
            datefmt="%H:%M:%S",
        )

    handler.setFormatter(formatter)
    root.addHandler(handler)

    # Keep noisy libraries sane
    logging.getLogger("werkzeug").setLevel(
        logging.WARNING if not debug else logging.INFO
    )
    logging.getLogger("dash").setLevel(logging.INFO if debug else logging.WARNING)

    # Make sure gunicorn logs propagate to root so they use the same handler
    for name in ("gunicorn", "gunicorn.error", "gunicorn.access"):
        lg = logging.getLogger(name)
        lg.handlers.clear()
        lg.propagate = True
        lg.setLevel(level)
