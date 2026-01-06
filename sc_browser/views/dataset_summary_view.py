from __future__ import annotations

from typing import Any, Dict

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from sc_browser.core.base_view import BaseView
from sc_browser.core.filter_state import FilterState
from sc_browser.core.filter_profile import FilterProfile


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
    filter_profile = FilterProfile(
        clusters=True,
        conditions=True,
        samples=True,
        cell_types=True,
        genes=False,
        embedding=False,
        split_by_condition=False,
        is_3d=False,
    )

    def compute_data(self, state: FilterState) -> Dict[str, Any]:

        # Apply filters via the Dataset abstraction (hits subset cache)
        ds = self.filtered_dataset(state)
        adata = ds.adata

        if adata.n_obs == 0:
            return {}

        n_cells, n_genes = adata.n_obs, adata.n_vars

        # obs schema table (column name, dtype, unique values)
        # This is mainly for debugging, so if it ever becomes expensive on huge datasets,
        # we could gate it behind a debug flag.
        obs_schema = adata.obs.dtypes.reset_index().rename(
            columns={"index": "column", 0: "dtype"}
        )
        obs_schema["n_unique"] = obs_schema["column"].map(
            lambda col: adata.obs[col].nunique()
        )

        # Cluster counts using pre-normalised Series from Dataset
        cluster_counts = ds.clusters.value_counts().reset_index()
        cluster_counts.columns = ["cluster", "count"]

        # Condition counts – only if a condition_key is configured
        if ds.condition_key is not None and ds.conditions is not None:
            condition_counts = ds.conditions.value_counts().reset_index()
            condition_counts.columns = ["condition", "count"]
        else:
            condition_counts = pd.DataFrame(columns=["condition", "count"])

        return {
            "n_cells": n_cells,
            "n_genes": n_genes,
            "obs_schema": obs_schema,
            "cluster_counts": cluster_counts,
            "condition_counts": condition_counts,
        }

    def render_figure(self, data: Dict[str, Any], state: FilterState) -> go.Figure:
        # Treat both None and {} as "no data"
        if not data:
            return self.empty_figure("No data to show")

        n_cells = data["n_cells"]
        n_genes = data["n_genes"]
        cluster_counts: pd.DataFrame = data["cluster_counts"]
        condition_counts: pd.DataFrame = data["condition_counts"]

        fig = make_subplots(
            rows=1,
            cols=2,
            subplot_titles=("Cluster sizes", "Condition sizes"),
        )

        if not cluster_counts.empty:
            fig.add_bar(
                x=cluster_counts["cluster"],
                y=cluster_counts["count"],
                row=1,
                col=1,
                name="Clusters",
            )

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
