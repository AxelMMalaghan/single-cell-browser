from __future__ import annotations

from typing import TYPE_CHECKING

import dash
from dash import Input, Output

from sc_browser.ui.ids import IDs
from sc_browser.ui.helpers import obs_preview_table

if TYPE_CHECKING:
    from sc_browser.ui.config import AppConfig


def register_dataset_preview_callbacks(app: dash.Dash, ctx: AppConfig) -> None:
    @app.callback(
        Output(IDs.Control.DM_OBS_PREVIEW, "children"),
        Input(IDs.Control.DATASET_SELECT, "value"),
    )
    def update_dataset_preview(dataset_name: str | None):
        # Dash can pass None/"" during initial load or when cleared
        if not dataset_name:
            return (
                "No dataset selected. "
                "Choose a dataset from the navbar dropdown to preview its .obs table."
            )

        ds = ctx.dataset_by_name.get(dataset_name)
        if ds is None:
            return (
                f"Dataset '{dataset_name}' not found. "
                "It may have been removed or failed to load. "
                "Choose another dataset."
            )

        return obs_preview_table(ds)
