from __future__ import annotations

from typing import Any

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from sc_browser.core.base_view import BaseView
from sc_browser.core.state import FilterState

class Dotplot(BaseView):
    """
    Dot plot:
    - x-axis: genes
    - y-axis: cluster (or cluster-condition)
    - dot-size %:of cells expressing the gene in that group
    - dot-color: relative mean expression in that group
    """

    id = "dotplot"
    label = "Dot Plot"

    def compute_data(self, state: FilterState) -> pd.DataFrame:
        # Apply filtering through Dataset abstraction
        ds = self.dataset.subset(
            clusters=state.clusters,
            conditions=state.conditions,
            samples=state.samples
        )
        adata_sub = ds.adata
        obs_sub = adata_sub.obs.copy()

        if adata_sub.n_obs == 0 or not state.genes:
            return pd.DataFrame()

        genes = [g for g in state.genes if g in adata_sub.var_names]
        if not genes:
            return pd.DataFrame()

        # --- expression matrix (cells x genes) ---
        expression_df = self.dataset.extract_expression_matrix(adata_sub, genes)
        if expression_df.empty:
            return pd.DataFrame()

        cluster_key = self.dataset.get_obs_column("cluster")
        condition_key = self.dataset.get_obs_column("condition")

        cluster = obs_sub[cluster_key].astype(str) if cluster_key else pd.Series("unknown", index=obs_sub.index)
        condition = obs_sub[condition_key].astype(str) if condition_key else pd.Series("unknown", index=obs_sub.index)

        if state.split_by_condition:
            cell_identity = cluster + "-" + condition
        else:
            cell_identity = cluster

        cell_meta = pd.DataFrame({
            "cellID": obs_sub.index,
            "cell_identity": cell_identity.values,
            "cluster": cluster.values,
            "condition": condition.values,
        })

        merged = long_expression.merge(cell_meta, on="cellID", how="left")

        def percent_expressed(x: pd.Series) -> float:
            return 100.0 * (x != 0).sum() / len(x) if len(x) else 0.0

        agg = (
            merged
            .groupby(["cell_identity", "gene"], as_index=False)["count"]
            .agg(meanExpr="mean", pctExpressed=percent_expressed)
        )

        agg["logMeanExpression"] = np.log2(1.0 + agg["meanExpr"])
        agg["maxMeanExpression"] = agg.groupby("gene")["meanExpr"].transform("max")
        agg["relMeanExpression"] = np.where(
            agg["maxMeanExpression"] > 0,
            agg["meanExpr"] / agg["maxMeanExpression"],
            0.0
        )

        return agg.sort_values(["gene", "cell_identity"]).reset_index(drop=True)

    def render_figure(self, data: pd.DataFrame, state: FilterState) -> Any:
        if data is None or data.empty:
            fig = go.Figure()
            fig.update_layout(
                title="Dot plot: select at least one gene to plot",
                xaxis={"visible": False},
                yaxis={"visible": False},
            )
            return fig

        fig = px.scatter(
            data,
            x="gene",
            y="cell_identity",
            size="pctExpressed",
            color="relMeanExpression",
            size_max=20,
            color_continuous_scale="viridis",
            hover_data={
                "pctExpressed": "mean",
                "meanExpr": "mean",
                "logMeanExpression": "mean",
            },
        )
        fig.update_layout(
            height=600,
            margin=dict(l=40, r=40, b=40, t=40),
            xaxis_title="Gene",
            yaxis_title="Cluster / Group",
            coloraxis_colorbar=dict(title = "Relative Mean Expression"),
        )
        fig.update_xaxes(tickangle=-45)

        return fig


