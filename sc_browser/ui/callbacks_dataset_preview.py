from __future__ import annotations

import dash

from dash import Input, Output
from typing import TYPE_CHECKING

from .helpers import obs_preview_table

if TYPE_CHECKING:
    from .config import AppConfig


def register_dataset_preview_callbacks(app: dash.Dash, ctx: AppConfig) -> None:
    @app.callback(
        Output("dm-obs-preview", "children"),
        Input("dataset-select", "value"),
    )
    def update_dataset_preview(dataset_name: str):
        if not dataset_name or dataset_name not in ctx.dataset_by_name:
            return (
                "No dataset selected. "
                "Choose a dataset from the navbar dropdown to preview its .obs table."
            )




        ds = ctx.dataset_by_name[dataset_name]
        return obs_preview_table(ds)
