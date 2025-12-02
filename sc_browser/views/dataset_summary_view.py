# sc_browser/views/dataset_summary_view.py

from __future__ import annotations

from typing import Any, Dict

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from sc_browser.core.base_view import BaseView
from sc_browser.core.state import FilterState
from sc_browser.core.dataset import Dataset


class DatasetSummary(BaseView):
    """
    Lightweight dataset inspector.

    Shows:
      - n_cells, n_genes
      - bar chart of cluster sizes
      - bar chart of condition sizes

    This is mainly a debugging / sanity-check view to confirm that:
      - config keys match AnnData.obs
      - filters (clusters/conditions) are being applied as expected
    """

    id = "dataset_summary"
    label = "Dataset Summary"

    def compute_data(self, state: FilterState) -> Dict[str, Any]:
        # Apply filters via the Dataset abstraction
        ds: Dataset = self.dataset.subset(
            clusters=state.clusters or None,
            conditions=state.conditions or None,
            samples=state.samples or None,
        )

        adata = ds.adata

        if adata.n_obs == 0:
            return {}



        n_cells, n_genes = adata.n_obs, adata.n_vars

        # obs schema table (column name, dtype, unique values)
        obs_schema = (
            adata.obs.dtypes.reset_index()
            .rename(columns={"index": "column", 0: "dtype"})
        )
        obs_schema["n_unique"] = obs_schema["column"].map(lambda col: adata.obs[col].nunique())

        # Cluster counts
        cluster_counts = ds.clusters.value_counts().reset_index()
        cluster_counts.columns = ["cluster", "count"]

        # Condition counts
        condition_counts = ds.conditions.value_counts().reset_index()
        condition_counts.columns = ["condition", "count"]

        return {
            "n_cells": n_cells,
            "n_genes": n_genes,
            "obs_schema": obs_schema,
            "cluster_counts": cluster_counts,
            "condition_counts": condition_counts,
        }

    def render_figure(self, data: Dict[str, Any], state: FilterState) -> go.Figure:
        # data is a dict, not a DataFrame – so do NOT call data.empty

        if not data:
            fig = go.Figure()
            fig.update_layout(
                title="No cells after filtering – adjust cluster/condition filters",
                xaxis={"visible": False},
                yaxis={"visible": False},
            )
            return fig

        n_cells = data["n_cells"]
        n_genes = data["n_genes"]
        cluster_counts: pd.DataFrame = data["cluster_counts"]
        condition_counts: pd.DataFrame = data["condition_counts"]

        # Build a 1x2 subplot: clusters (left), conditions (right)
        fig = make_subplots(
            rows=1,
            cols=2,
            subplot_titles=("Cluster sizes", "Condition sizes"),
        )

        # Cluster bar chart
        if not cluster_counts.empty:
            fig.add_bar(
                x=cluster_counts["cluster"],
                y=cluster_counts["count"],
                row=1,
                col=1,
                name="Clusters",
            )

        # Condition bar chart
        if not condition_counts.empty:
            fig.add_bar(
                x=condition_counts["condition"],
                y=condition_counts["count"],
                row=1,
                col=2,
                name="Conditions",
            )

        fig.update_xaxes(title_text="Cluster", row=1, col=1)
        fig.update_yaxes(title_text="# cells", row=1, col=1)

        fig.update_xaxes(title_text="Condition", row=1, col=2)
        fig.update_yaxes(title_text="# cells", row=1, col=2)

        fig.update_layout(
            height=600,
            margin=dict(l=40, r=40, t=60, b=40),
            title=f"Dataset summary: {n_cells} cells, {n_genes} genes",
            showlegend=False,
        )

        return fig
