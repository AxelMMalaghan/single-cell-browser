from __future__ import annotations

from typing import Any, Dict, List, Sequence

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from sc_browser.core.dataset import Dataset
from sc_browser.core.filter_state import FilterState
from sc_browser.core.filter_profile import FilterProfile
from sc_browser.core.base_view import BaseView


class DotplotView(BaseView):
    """
    Dot plot of gene expression across groups.

    - x: gene
    - y: cell_identity (cluster or cluster|condition)
    - dot size: pctExpressed (percent of cells expressing gene in group)
    - dot color: relMeanExpression (mean expression scaled per gene)
    """

    id = "dotplot"
    label = "Dot plot"
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

    def _cell_identities(self, ds: Dataset, state: FilterState) -> pd.Series:
        """
        Compute the 'cell_identity' label per cell, based on cluster and (optional) condition.
        """
        clusters = ds.clusters  # pre-normalised in Dataset
        conditions = ds.conditions

        # If split_by_condition is enabled and conditions exist, combine them
        if getattr(state, "split_by_condition", False) and conditions is not None:
            return clusters.astype(str) + " | " + conditions.astype(str)

        # Otherwise, just use cluster as the identity
        return clusters.astype(str)

    def compute_data(self, state: FilterState) -> pd.DataFrame:
        ds: Dataset = self.dataset

        genes: List[str] = state.genes or []
        # No genes selected so nothing to plot
        if not genes:
            return pd.DataFrame()

        # Apply filters (uses subset cache)
        ds_sub = self.filtered_dataset(state)

        adata = ds_sub.adata
        if adata.n_obs == 0:
            return pd.DataFrame()

        # Expression matrix for selected genes (uses expression cache)
        expr = ds_sub.expression_matrix(genes)  # cells × genes, dense

        if expr.empty:
            return pd.DataFrame()

        # Build long-form table: cell, gene, expr
        expr_long = expr.stack().rename("expr").reset_index()
        expr_long.columns = ["cell", "gene", "expr"]

        # Attach group labels
        expr_long["cluster"] = ds_sub.clusters.loc[expr_long["cell"]].values
        cell_identity = self._cell_identities(ds_sub, state)
        expr_long["cell_identity"] = cell_identity.loc[expr_long["cell"]].values

        if ds_sub.condition_key is not None and ds_sub.conditions is not None:
            expr_long["condition"] = ds_sub.conditions.loc[expr_long["cell"]].values
        else:
            expr_long["condition"] = None

        # Group by (cell_identity, gene, cluster, condition)
        group_cols = ["cell_identity", "gene", "cluster", "condition"]

        grouped = expr_long.groupby(group_cols, as_index=False).agg(
            n_cells=("expr", "size"),
            n_expressing=("expr", lambda x: (x > 0).sum()),
            meanExpr=("expr", "mean"),
        )

        # Percent expressing
        grouped["pctExpressed"] = (
            100.0
            * grouped["n_expressing"]
            / grouped["n_cells"].where(grouped["n_cells"] != 0, np.nan)
        )

        # Relative mean expression per gene (scale within each gene)
        # relMeanExpression = meanExpr / mean(meanExpr for that gene)
        gene_mean = grouped.groupby("gene")["meanExpr"].transform("mean")
        grouped["relMeanExpression"] = grouped["meanExpr"] / (
            gene_mean.replace(0, np.nan)
        )

        # Log mean expression for hover info
        grouped["logMeanExpression"] = np.log1p(grouped["meanExpr"].clip(lower=0))

        # Keep only the columns used by render_figure
        cols = [
            "gene",
            "cell_identity",
            "pctExpressed",
            "relMeanExpression",
            "meanExpr",
            "logMeanExpression",
            "cluster",
            "condition",
        ]
        return grouped[cols]

    def render_figure(self, data: pd.DataFrame, state: FilterState) -> Any:

        if data is None or data.empty:
            return self.empty_figure("No data to show")

        color_scale = state.color_scale

        fig = px.scatter(
            data,
            x="gene",
            y="cell_identity",
            size="pctExpressed",
            color="relMeanExpression",
            size_max=20,
            color_continuous_scale=color_scale,
            hover_data={
                "pctExpressed": True,
                "meanExpr": True,
                "logMeanExpression": True,
                "cluster": True,
                "condition": True,
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
