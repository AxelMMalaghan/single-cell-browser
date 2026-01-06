from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any, Optional

import dash
import pandas as pd
import plotly.graph_objs as go
from dash import Input, Output

from sc_browser.core.filter_state import FilterState
from sc_browser.ui.ids import IDs

if TYPE_CHECKING:
    from sc_browser.ui.config import AppConfig

logger = logging.getLogger(__name__)


# -----------------------------------------------------------------------------
# Helper: Empty/Error Figures
# -----------------------------------------------------------------------------
def _message_figure(title: str, details: Optional[str] = None) -> go.Figure:
    fig = go.Figure()
    text = title if details is None else f"{title}\n\n{details}"
    fig.add_annotation(
        text=text,
        showarrow=False,
        xref="paper",
        yref="paper",
        x=0.5,
        y=0.5,
    )
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    fig.update_layout(margin=dict(l=40, r=40, t=40, b=40))
    return fig


def _error_figure(details: str) -> go.Figure:
    return _message_figure("Something went wrong while rendering this view.", details)


def register_render_callbacks(app: dash.Dash, ctx: AppConfig) -> None:
    # ---------------------------------------------------------
    # Main figure: FilterState -> figure
    # ---------------------------------------------------------
    @app.callback(
        Output(IDs.Control.MAIN_GRAPH, "figure"),
        Input(IDs.Store.FILTER_STATE, "data"),
        prevent_initial_call=True,
    )
    def update_main_graph_from_state(fs_data: dict[str, Any] | None):
        if fs_data is None:
            return _message_figure(
                "No dataset/view selected.",
                "Choose a dataset and view to see any plots.",
            )

        try:
            state = FilterState.from_dict(fs_data)
        except Exception:
            logger.exception("Invalid filter state in main graph callback: %r", fs_data)
            return _error_figure("Internal error: invalid filter state.")

        ds = ctx.dataset_by_name.get(state.dataset_name)
        if ds is None:
            return _error_figure(
                f"The dataset '{state.dataset_name}' is not available. "
                "Try reloading the app or selecting a different dataset."
            )

        # Basic pre-flight checks (optional optimization to avoid loading data)
        if fs_data.get("genes") and not state.genes:
            return _message_figure(
                "Selected genes not found in this dataset.",
                "Try a different gene symbol or clear the gene filter.",
            )

        if state.embedding and (state.dim_x is None or state.dim_y is None):
            return _message_figure(
                "Select dimensions for the embedding.",
                "Choose X and Y dimensions (and Z for 3D, if enabled) from the sidebar.",
            )

        try:
            # 1. Instantiate View
            registry = ctx.registry
            if registry is None:
                return _error_figure("View registry is not available.")
            view = registry.create(state.view_id, ds)

            logger.info(
                "render_start",
                extra={
                    "view_id": state.view_id,
                    "dataset": state.dataset_name,
                    "n_genes": len(state.genes),
                },
            )

            # 2. Compute Data (Timed)
            data = view.timed_compute(state)

            if isinstance(data, pd.DataFrame) and data.empty:
                return _message_figure(
                    "No data to display.",
                    "Your current filters removed all cells for this view. "
                    "Try clearing one or more filters (clusters, conditions, samples, or cell types).",
                )

            if data is None:
                return _message_figure(
                    "No data returned by this view.",
                    "This can happen if there are no matching cells or the view "
                    "does not support the current settings.",
                )

            # 3. Render Figure
            return view.render_figure(data, state)

        except Exception:
            logger.exception(
                "Error in update_main_graph_from_state",
                extra={"filter_state": fs_data},
            )
            return _error_figure(
                "The app hit an unexpected error. "
                "If this keeps happening, grab the logs and open an issue."
            )
