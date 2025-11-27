from __future__ import annotations

from typing import Any

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from sc_browser.core.view_base import BaseView
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
        adata = self.dataset.adata
        obs = adata.obs.copy()

        # --- mask by cluster / condition
        mask = np.ones(len(obs), dtype=bool)

        if state.clusters:
            mask &= obs[self.dataset.cluster_key].astype(str).isin(state.clusters)

        if state.conditions:
            mask &= obs[self.dataset.condition_key].astype(str).isin(state.conditions)

        adata_sub = adata[mask].copy()
        obs_sub = obs[mask].copy()

        if adata_sub.n_obs == 0:
            return pd.DataFrame()

        # --- gene selection ---

        if not state.genes or len(state.genes) == 0:
            return pd.DataFrame()

        genes = [g for g in state.genes if g in adata_sub.var_names]
        if len(genes) == 0:
            return pd.DataFrame()

        # --- matrix (cell x genes) ---
        X = adata_sub[:, genes].X
        if hasattr(X, "toarray"):
            X = X.toarray()

        expression_df = pd.DataFrame(
            X,
            index=adata_sub.obs_names,
            columns=genes
        )

        # long format: cellID, gene, count
        long_expression = (
            expression_df
            .reset_index(names="cellID")
            .melt(id_vars="cellID", var_name="gene", value_name="count")
        )

        # --- group (cluster or cluster-condition) ---
        cluster = obs_sub[self.dataset.cluster_key].astype(str)
        condition = obs_sub[self.dataset.condition_key].astype(str)

        if state.split_by_condition:
            cell_identity = cluster + "-" + condition
        else:
            cell_identity = cluster

        cell_meta = pd.DataFrame(
            {
                "cellID": obs_sub.index,
                "cell_identity": cell_identity.values,
                "cluster": cluster.values,
                "condition": condition.values,
            }
        )

        merged = long_expression.merge(cell_meta, on="cellID", how="left")

        # --- aggregate per cell ---
        def percent_expressed(x: pd.Series) -> float:
            if x.sum() == 0:
                return np.nan
            return 100.0 * (x != 0).sum() / len(x)

        grouped = merged.groupby(["cell_identity", "gene"], as_index=False)
        agg = grouped["count"].agg(
            meanExpr="mean",
            pctExpressed=percent_expressed,
        )

        # log mean expression
        agg["logMeanExpression"] = np.log2(1.0 + agg["meanExpr"])

        #relative expression per gene
        agg["maxMeanExpression"] = agg.groupby("gene")["meanExpr"].transform(
            lambda s: s.max(skipna=True)
        )
        agg["relMeanExpression"] = np.where(
            agg["maxMeanExpression"] > 0,
            agg["meanExpr"] / agg["maxMeanExpression"],
            0.0
        )

        agg = agg.sort_values(["gene", "cell_identity"]).reset_index(drop=True)

        return agg

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


