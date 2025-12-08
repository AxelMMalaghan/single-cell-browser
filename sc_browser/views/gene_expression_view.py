from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go

from sc_browser.core.filter_state import FilterState, FilterProfile
from sc_browser.core.base_view import BaseView
from sc_browser.core.dataset import Dataset


class ExpressionView(BaseView):
    """
    Expression plot:
    Embed cells in 2D and colour by log1p(gene expression).
    """

    id = "expression"
    label = "Expression"
    filter_profile = FilterProfile(
        clusters=True,
        conditions=True,
        samples=True,
        cell_types=True,
        genes=True,
        embedding=False
    )

    def compute_data(self, state: FilterState) -> pd.DataFrame:
        ds: Dataset = self.dataset

        # --- gene selection ---
        if not state.genes:
            return pd.DataFrame()

        # MVP: only the first gene
        gene = state.genes[0]

        # Apply filters using Dataset abstraction (hits subset cache)
        ds_sub = self.filtered_dataset(state)

        adata = ds_sub.adata
        if adata.n_obs == 0:
            return pd.DataFrame()

        if gene not in adata.var_names:
            return pd.DataFrame()

        # Use Dataset.expression_matrix for dense extraction + caching
        expr_df = ds_sub.expression_matrix([gene])
        if expr_df.empty:
            return pd.DataFrame()

        expression = expr_df[gene].to_numpy()
        obs_index = expr_df.index

        # --- embedding ---
        emb_key = state.embedding or ds.embedding_key
        embedding = ds_sub.get_embedding(emb_key).loc[obs_index].copy()
        embedding = embedding.rename(columns={"dim1": "x", "dim2": "y"})

        # --- grouping ---
        clusters = ds_sub.clusters.loc[obs_index] if ds_sub.clusters is not None else pd.Series(
            ["NA"] * len(obs_index),
            index=obs_index,
        )

        conditions = (
            ds_sub.conditions.loc[obs_index]
            if ds_sub.condition_key is not None and ds_sub.conditions is not None
            else pd.Series(["NA"] * len(obs_index), index=obs_index)
        )

        if state.split_by_condition:
            group = clusters.astype(str) + "_" + conditions.astype(str)
        else:
            group = clusters.astype(str)

        df = pd.DataFrame(
            {
                "x": embedding["x"],
                "y": embedding["y"],
                "expression": expression,
                "log_expression": np.log1p(expression),
                "cluster": clusters.astype(str),
                "condition": conditions.astype(str),
                "group": group.astype(str),
                "gene": gene,
            },
            index=obs_index,
        )

        return df

    def render_figure(self, data: pd.DataFrame, state: FilterState) -> go.Figure:

        if data is None or data.empty:
            return self.empty_figure("No data to show")


        # Data plotting
        color_scale = getattr(state, "color_scale", "viridis")

        fig = px.scatter(
            data,
            x="x",
            y="y",
            color="log_expression",
            color_continuous_scale=color_scale,
            hover_data={
                "cluster": True,
                "condition": True,
                "expression": True,
                "log_expression": True,
                "gene": True,
            },
        )

        fig.update_layout(
            height=600,
            margin=dict(l=40, r=40, b=40, t=40),
            coloraxis_colorbar=dict(title="log1p(expression)"),
        )

        return fig
