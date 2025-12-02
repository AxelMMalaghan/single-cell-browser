from __future__ import annotations

from typing import Any, List

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from sc_browser.core.state import FilterState
from sc_browser.core.base_view import BaseView

class HeatmapView(BaseView):
    """
    Heatmap of average expression per group (cluster or cluster-condition)
    for the selected genes
    """

    id = "heatmap"
    label = "Heatmap"

    def compute_data(self, state: FilterState) -> pd.DataFrame:
        # Apply filters using the unified Dataset.subset abstraction
        ds = self.dataset.subset(
            clusters=state.clusters or None,
            conditions=state.conditions or None,
            samples=state.samples or None,
        )

        adata = ds.adata

        if adata.n_obs == 0 or not state.genes:
            return pd.DataFrame()

        genes = [g for g in state.genes if g in adata.var_names]
        if not genes:
            return pd.DataFrame()

            # Filter genes to those present in the dataset
        genes: List[str] = [g for g in state.genes if g in adata.var_names]
        if not genes:
            return pd.DataFrame()

        # cells x genes matrix
        expression_df = self.dataset.extract_expression_matrix(adata, genes)
        if expression_df.empty:
            return pd.DataFrame()

        cluster = ds.clusters.astype(str)
        condition = ds.conditions.astype(str)

        group = cluster + "_" + condition if state.split_by_condition else cluster
        expression_df["group"] = group.values

        # Compute average expression per group and gene
        long_df = (
            expression_df
            .melt(id_vars="group", var_name="gene", value_name="expression")
            .groupby(["group", "gene"], as_index=False)
            .agg(mean_expression=("expression", "mean"))
        )

        long_df["log_mean_expression"] = np.log1p(long_df["mean_expression"])

        long_df["gene"] = pd.Categorical(
            long_df["gene"],
            categories=[g for g in genes if g in long_df["gene"].unique()],
            ordered=True,
        )

        long_df["group"] = pd.Categorical(
            long_df["group"],
            categories=sorted(long_df["group"].unique()),
            ordered=True,
        )

        return long_df

    def render_figure(self, data: pd.DataFrame, state: FilterState) -> Any:

        if data is None or data.empty:
            fig = go.Figure()
            fig.update_layout(
                title="Select at least one valid gene to show heatmap",
                xaxis={"visible": False},
                yaxis={"visible": False},
            )
            return fig

        # pivot to matrix genes x genes
        pivot = data.pivot(index="gene", columns="group", values="log_mean_expression")

        fig = px.imshow(
            pivot,
            color_continuous_scale="Viridis",
            aspect="auto",
            labels=dict(x="Group", y="Gene", color="log1p(mean expression)"),
        )

        fig.update_xaxes(side="top")
        fig.update_layout(
            height=600,
            margin=dict(l=40, r=40, b=40, t=40),
        )

        return fig

