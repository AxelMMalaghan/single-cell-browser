from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go

from sc_browser.core.state import FilterState
from sc_browser.core.base_view import BaseView


class ExpressionView(BaseView):
    """
    Expression plot:
    Embed cells in 2D and colour by log1p(gene expression).
    """

    id = "expression"
    label = "Expression"

    def compute_data(self, state: FilterState) -> pd.DataFrame:
        ds = self.dataset.subset(
            clusters=state.clusters or None,
            conditions=state.conditions or None,
            samples=state.samples or None,
        )

        adata = ds.adata

        # --- gene selection ---
        if not state.genes or len(state.genes) == 0:
            return pd.DataFrame()

        gene = state.genes[0]  # MVP: only first gene

        if gene not in adata.var_names:
            return pd.DataFrame()

        expression_df = ds.extract_expression_matrix(adata, [gene])
        if expression_df.empty:
            return pd.DataFrame()

        expression = expression_df[gene].to_numpy()
        obs_index = expression_df.index

        # --- embedding ---
        embedding = ds.embedding.loc[obs_index].copy()
        embedding = embedding.rename(columns={"dim1": "x", "dim2": "y"})

        # --- grouping ---
        cluster = ds.clusters.loc[obs_index].astype(str)
        condition = ds.conditions.loc[obs_index].astype(str)

        group = cluster + "_" + condition if state.split_by_condition else cluster

        df = pd.DataFrame({
            "x": embedding["x"],
            "y": embedding["y"],
            "expression": expression,
            "cluster": cluster,
            "condition": condition,
            "group": group,
            "gene": gene,
        }, index=obs_index)

        df["log_expression"] = np.log1p(df["expression"])
        return df

    def render_figure(self, data: pd.DataFrame, state: FilterState) -> Any:
        if data.empty:
            fig = go.Figure()
            fig.update_layout(
                title="Select at least one gene to show expression",
                xaxis={"visible": False},
                yaxis={"visible": False},
            )
            return fig

        fig = px.scatter(
            data,
            x="x",
            y="y",
            color="log_expression",
            color_continuous_scale="viridis",
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
