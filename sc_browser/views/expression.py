from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
from sc_browser.core.state import FilterState
from sc_browser.core.view_base import BaseView

class ExpressionView(BaseView):
    """
    Expression plot:
    embed cells in 2d and colour by log1p (gene expression)
    """

    id = "expression"
    label = "Expression"

    def compute_data(self, state: FilterState) -> pd.DataFrame:

        adata = self.dataset.adata
        obs = adata.obs.copy()

        # --- masking by cluster / condition ---
        mask = np.ones(len(obs), dtype=bool)

        if state.clusters:
            mask &= obs[self.dataset.cluster_key].astype(str).isin(state.clusters)

        if state.conditions:
            mask &= obs[self.dataset.condition_key].astype(str).isin(state.conditions)

        adata_sub = adata[mask].copy()
        obs_sub = adata_sub.obs.copy()

        # --- gene selection ---
        if not state.genes or len(state.genes) == 0:
            return pd.DataFrame()

        gene = state.genes[0]   #MVP: first gene only

        if gene not in adata_sub.var_names:
            #Unknown gene -> empty
            return pd.DataFrame()

        X = adata_sub[:, gene].X
        if hasattr(X, "toarray"):
            X = X.toarray()
        expression = np.asarray(X).ravel()


        # --- embedding ---
        embedding = self.dataset.embedding.loc[obs_sub.index].copy()
        embedding = embedding.rename(columns={"dim1": "x", "dim2": "y"})

        # --- grouping / facet concept ----
        cluster = obs_sub[self.dataset.cluster_key].astype(str)
        condition = obs_sub[self.dataset.condition_key].astype(str)


        if state.split_by_condition:
            group = cluster + "_" + condition
        else:
            group = cluster

        df = pd.DataFrame(
            {
                "x": embedding["x"],
                "y": embedding["y"],
                "expression": expression,
                "cluster": cluster,
                "condition": condition,
                "group": group,
                "gene": gene,
            },
            index=obs_sub.index,
        )

        # Precompute log expression to send to figure
        df["log_expression"] = np.log1p(df["expression"])

        return df

    def render_figure(self, data: pd.DataFrame, state: FilterState) -> Any:

        if data.empty:
            # Return empty figure
            fig = go.Figure()
            fig.update_layout(
                title="Select atleast one gene to show expression",
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

