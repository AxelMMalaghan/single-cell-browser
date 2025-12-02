from __future__ import annotations

from typing import Any

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from sc_browser.core.state import FilterState
from sc_browser.core.base_view import BaseView


class FeatureCountView(BaseView):
    """
    QC plot: per cell total counts vs number of detected features
    """

    id = "feature_count"
    label = "Feature v. Count"

    def compute_data(self, state: FilterState) -> pd.DataFrame:
        # Apply filters using Dataset abstraction
        ds = self.dataset.subset(
            clusters=state.clusters or None,
            conditions=state.conditions or None,
            samples=state.samples or None,
        )

        adata = ds.adata


        if adata.n_obs == 0:
            return pd.DataFrame()

        genes = list(adata.var_names)

        # Use unified method to get expression matrix
        expression_df = ds.extract_expression_matrix(adata, genes)
        if expression_df.empty:
            return pd.DataFrame()

        n_counts = expression_df.sum(axis=1)
        n_features = (expression_df > 0).sum(axis=1)

        cluster = ds.clusters.astype(str)
        condition = ds.conditions.astype(str)

        group = cluster + "_" + condition if state.split_by_condition else cluster

        df = pd.DataFrame({
            "n_counts": n_counts.values,
            "n_features": n_features.values,
            "cluster": cluster.values,
            "condition": condition.values,
            "group": group.values,
        }, index=expression_df.index)

        return df

    def render_figure(self, data: pd.DataFrame, state: FilterState) -> Any:
        if data is None or data.empty:
            fig = go.Figure()
            fig.update_layout(
                title="No cells after filtering - adjust cluster/condition filters",
                xaxis={"visible": False},
                yaxis={"visible": False},
            )
            return fig

        fig = px.scatter(
            data,
            x="n_counts",
            y="n_features",
            color="group",
            hover_data={"cluster": True, "condition": True, "n_counts": True, "n_features": True},
        )

        fig.update_layout(
            height=600,
            margin=dict(l=40, r=40, t=40, b=40),
            xaxis_title="Total counts per cell",
            yaxis_title="Number of features per cell",
            legend_title="Group"
        )
        return fig
