"""
UI adapters for the browser.

Currently provides a Dash-based web UI via create_dash_app().
Other adapters (REST API, Jupyter, etc.) can live here later.
"""

from .dash_app import create_dash_app

__all__ = ["create_dash_app"]
