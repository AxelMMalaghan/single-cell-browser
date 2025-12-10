from __future__ import annotations

from typing import List

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from sc_browser.core.filter_state import FilterState
from sc_browser.core.filter_profile import FilterProfile
from sc_browser.core.base_view import BaseView
from sc_browser.core.dataset import Dataset


class HeatmapView(BaseView):
    """
    Heatmap of average expression per group (cluster or cluster-condition)
    for the selected genes.
    """

    id = "heatmap"
    label = "Heatmap"
    filter_profile = FilterProfile(
        clusters=True,
        conditions=True,
        samples=True,
        cell_types=True,
        genes=True,
        embedding=False,
        split_by_condition=True,
        is_3d=False,
        colour_scale=True,
    )

    def compute_data(self, state: FilterState) -> pd.DataFrame:
        base_ds: Dataset = self.dataset

        # No genes selected → nothing to do
        if not state.genes:
            return pd.DataFrame()

        # Apply filters using the unified Dataset.subset abstraction (hits subset cache)
        ds = self.filtered_dataset(state)

        adata = ds.adata
        if adata.n_obs == 0:
            return pd.DataFrame()

        # Filter genes to those present in the dataset
        genes: List[str] = [g for g in state.genes if g in adata.var_names]
        if not genes:
            return pd.DataFrame()

        # cells × genes expression matrix (dense, cached)
        expr_df = ds.expression_matrix(genes)
        if expr_df.empty:
            return pd.DataFrame()

        # Pre-normalised labels from Dataset
        clusters = ds.clusters
        if clusters is None:
            clusters = pd.Series(["NA"] * adata.n_obs, index=adata.obs_names)

        conditions = (
            ds.conditions
            if ds.condition_key is not None and ds.conditions is not None
            else pd.Series(["NA"] * adata.n_obs, index=adata.obs_names)
        )

        if state.split_by_condition:
            group = clusters.astype(str) + "_" + conditions.astype(str)
        else:
            group = clusters.astype(str)

        # Attach group to expression_df
        expr_df = expr_df.copy()
        expr_df["group"] = group.loc[expr_df.index].values

        # Long-form: group, gene, expression
        long_df = (
            expr_df
            .melt(id_vars="group", var_name="gene", value_name="expression")
            .groupby(["group", "gene"], as_index=False)
            .agg(mean_expression=("expression", "mean"))
        )

        # log1p(mean expression)
        long_df["log_mean_expression"] = np.log1p(long_df["mean_expression"])

        # Keep gene order as in the input list, but only those that survived
        ordered_genes = [g for g in genes if g in long_df["gene"].unique()]
        long_df["gene"] = pd.Categorical(
            long_df["gene"],
            categories=ordered_genes,
            ordered=True,
        )

        # Order groups alphabetically (or keep as-is if you prefer)
        long_df["group"] = pd.Categorical(
            long_df["group"],
            categories=sorted(long_df["group"].unique()),
            ordered=True,
        )

        return long_df

    def render_figure(self, data: pd.DataFrame, state: FilterState) -> go.Figure:

        if data is None or data.empty:
            return self.empty_figure("No data to show")

        # Pivot to matrix: genes × groups
        pivot = data.pivot(index="gene", columns="group", values="log_mean_expression")

        color_scale = state.color_scale

        fig = px.imshow(
            pivot,
            color_continuous_scale=color_scale,
            aspect="auto",
            labels=dict(x="Group", y="Gene", color="log1p(mean expression)"),
        )

        fig.update_xaxes(side="top")
        fig.update_layout(
            height=600,
            margin=dict(l=40, r=40, b=40, t=40),
        )

        return fig
