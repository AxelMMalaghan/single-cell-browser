from __future__ import annotations

from typing import Any

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from sc_browser.core.base_view import BaseView
from sc_browser.core.state import FilterState


class DotplotView(BaseView):
    """
    Dot plot:
    - x-axis: genes
    - y-axis: cluster (or cluster-condition)
    - dot-size: % of cells expressing the gene in that group
    - dot-color: relative mean expression in that group
    """

    id = "dotplot"
    label = "Dot Plot"

    def compute_data(self, state: FilterState) -> pd.DataFrame:
        # ---- Subset dataset according to filters ----
        samples = getattr(state, "samples", None)  # be robust if samples isn't wired yet

        ds = self.dataset.subset(
            clusters=state.clusters,
            conditions=state.conditions,
            samples=samples,
        )

        adata_sub = ds.adata
        obs_sub = adata_sub.obs.copy()

        # No cells or no genes selected
        if adata_sub.n_obs == 0 or not state.genes:
            return pd.DataFrame()

        # Keep only genes that exist in this dataset
        genes = [g for g in state.genes if g in adata_sub.var_names]
        if not genes:
            return pd.DataFrame()

        # --- expression matrix (cells x genes) ---
        # Assuming this returns a DataFrame with index = cell IDs, columns = gene names
        expression_df = ds.extract_expression_matrix(adata_sub, genes)
        if expression_df.empty:
            return pd.DataFrame()

        # Ensure we have a proper cellID index name for merging later
        expression_df = expression_df.copy()
        if expression_df.index.name != "cellID":
            expression_df.index = expression_df.index.astype(str)
            expression_df.index.name = "cellID"

        # Wide -> long: one row per (cell, gene)
        long_expression = (
            expression_df
            .reset_index()  # brings cellID out as a column
            .melt(
                id_vars="cellID",
                var_name="gene",
                value_name="count",
            )
        )

        # --- group definitions (cluster / condition) ---
        cluster_key = ds.get_obs_column("cluster")
        condition_key = ds.get_obs_column("condition")

        if cluster_key and cluster_key in obs_sub.columns:
            cluster = obs_sub[cluster_key].astype(str)
        else:
            cluster = pd.Series("unknown", index=obs_sub.index)

        if condition_key and condition_key in obs_sub.columns:
            condition = obs_sub[condition_key].astype(str)
        else:
            condition = pd.Series("unknown", index=obs_sub.index)

        if state.split_by_condition:
            cell_identity = cluster + "-" + condition
        else:
            cell_identity = cluster

        cell_meta = pd.DataFrame(
            {
                "cellID": obs_sub.index.astype(str),
                "cell_identity": cell_identity.values,
                "cluster": cluster.values,
                "condition": condition.values,
            }
        )

        # --- join expression with metadata ---
        merged = long_expression.merge(cell_meta, on="cellID", how="left")

        # --- aggregate per (group, gene) ---
        def percent_expressed(x: pd.Series) -> float:
            return 100.0 * (x != 0).sum() / len(x) if len(x) else 0.0

        agg = (
            merged
            .groupby(["cell_identity", "gene"], as_index=False)["count"]
            .agg(meanExpr="mean", pctExpressed=percent_expressed)
        )

        # --- normalise expression for coloring ---
        agg["logMeanExpression"] = np.log2(1.0 + agg["meanExpr"])
        agg["maxMeanExpression"] = agg.groupby("gene")["meanExpr"].transform("max")
        agg["relMeanExpression"] = np.where(
            agg["maxMeanExpression"] > 0,
            agg["meanExpr"] / agg["maxMeanExpression"],
            0.0,
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
                "pctExpressed": True,
                "meanExpr": True,
                "logMeanExpression": True,
            },
        )

        fig.update_layout(
            height=600,
            margin=dict(l=40, r=40, b=40, t=40),
            xaxis_title="Gene",
            yaxis_title="Cluster / Group",
            coloraxis_colorbar=dict(title="Relative Mean Expression"),
        )
        fig.update_xaxes(tickangle=-45)

        return fig