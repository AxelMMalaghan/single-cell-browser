from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from sc_browser.core.base_view import BaseView
from sc_browser.core.filter_state import FilterState, FilterProfile
from sc_browser.core.dataset import Dataset


class FeatureCountView(BaseView):
    """
    QC plot: per cell total counts vs number of detected features.
    """

    id = "feature_count"
    label = "Feature v. Count"
    filter_profile = FilterProfile()

    def compute_data(self, state: FilterState) -> pd.DataFrame:
        # Apply filters using Dataset abstraction (hits subset cache)
        ds: Dataset = self.dataset.subset(
            clusters=state.clusters or None,
            conditions=state.conditions or None,
            samples=getattr(state, "samples", None) or None,
            cell_types=getattr(state, "cell_types", None) or None,
        )

        adata = ds.adata

        if adata.n_obs == 0:
            return pd.DataFrame()

        # Work directly on the matrix: cells x genes
        X = adata.X

        # Total counts per cell (works for dense and sparse)
        n_counts = np.asarray(X.sum(axis=1)).ravel()

        # Number of detected features per cell: how many genes > 0
        n_features = np.asarray((X > 0).sum(axis=1)).ravel()

        # Pre-normalised labels from Dataset
        cluster = ds.clusters
        if cluster is None:
            cluster = pd.Series(["NA"] * adata.n_obs, index=adata.obs_names)

        conditions = ds.conditions
        if conditions is None:
            conditions = pd.Series(["NA"] * adata.n_obs, index=adata.obs_names)

        # Group for coloring: cluster or cluster_condition if split_by_condition
        if state.split_by_condition:
            group = cluster.astype(str) + "_" + conditions.astype(str)
        else:
            group = cluster.astype(str)

        df = pd.DataFrame(
            {
                "n_counts": n_counts,
                "n_features": n_features,
                "cluster": cluster.values,
                "condition": conditions.values,
                "group": group.values,
            },
            index=adata.obs_names,
        )

        return df

    def render_figure(self, data: pd.DataFrame, state: FilterState) -> Any:
        if data is None or data.empty:
            fig = go.Figure()
            fig.update_layout(
                title="No cells after filtering - adjust filters",
                xaxis={"visible": False},
                yaxis={"visible": False},
            )
            return fig

        fig = px.scatter(
            data,
            x="n_counts",
            y="n_features",
            color="group",
            hover_data={
                "cluster": True,
                "condition": True,
                "n_counts": True,
                "n_features": True,
            },
        )

        fig.update_layout(
            height=600,
            margin=dict(l=40, r=40, t=40, b=40),
            xaxis_title="Total counts per cell",
            yaxis_title="Number of features per cell",
            legend_title="Group",
        )
        return fig
