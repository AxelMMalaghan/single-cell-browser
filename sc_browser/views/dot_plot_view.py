from __future__ import annotations

from typing import Any, List

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scipy.sparse as sp

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
        samples = getattr(state, "samples", None)

        ds = self.dataset.subset(
            clusters=state.clusters or None,
            conditions=state.conditions or None,
            samples=samples or None,
        )

        adata = ds.adata

        # No cells or no genes selected
        if adata.n_obs == 0 or not state.genes:
            return pd.DataFrame()

        # Keep only genes that exist in this dataset
        genes: List[str] = [g for g in state.genes if g in adata.var_names]
        if not genes:
            return pd.DataFrame()

        # --- expression matrix (cells x genes) ---
        X = adata[:, genes].X
        if sp.issparse(X):
            X = X.tocsr()

        # --- group definitions (cluster / condition) ---
        obs = adata.obs

        cluster_key = ds.get_obs_column("cluster")
        condition_key = ds.get_obs_column("condition")

        if cluster_key and cluster_key in obs.columns:
            cluster = obs[cluster_key].astype(str)
        else:
            cluster = pd.Series("unknown", index=obs.index)

        if condition_key and condition_key in obs.columns:
            condition = obs[condition_key].astype(str)
        else:
            condition = pd.Series("unknown", index=obs.index)

        if state.split_by_condition:
            cell_identity = cluster + "-" + condition
        else:
            cell_identity = cluster

        # Categorical groups so we can index them efficiently
        groups = pd.Categorical(cell_identity)
        group_codes = groups.codes  # int code per cell
        group_names = list(groups.categories)
        n_groups = len(group_names)
        n_genes = len(genes)

        # Pre-allocate matrices: group x gene
        mean_expr = np.zeros((n_groups, n_genes), dtype=float)
        pct_expr = np.zeros((n_groups, n_genes), dtype=float)

        for g_idx in range(n_groups):
            mask = (group_codes == g_idx)
            if not np.any(mask):
                continue

            Xg = X[mask]

            if sp.issparse(Xg):
                # sums and non-zero counts per gene
                sums = np.asarray(Xg.sum(axis=0)).ravel()
                nonzero = np.asarray((Xg > 0).sum(axis=0)).ravel()
            else:
                sums = Xg.sum(axis=0)
                nonzero = (Xg > 0).sum(axis=0)

            n_cells = mask.sum()
            mean_expr[g_idx, :] = sums / n_cells
            pct_expr[g_idx, :] = (nonzero / n_cells) * 100.0

        # Build tidy DataFrame: one row per (group, gene)
        cell_identity_col = np.repeat(group_names, n_genes)
        gene_col = np.tile(genes, n_groups)

        df = pd.DataFrame(
            {
                "cell_identity": cell_identity_col,
                "gene": gene_col,
                "meanExpr": mean_expr.ravel(),
                "pctExpressed": pct_expr.ravel(),
            }
        )

        # Add cluster/condition for hover
        meta = (
            pd.DataFrame(
                {
                    "cell_identity": cell_identity,
                    "cluster": cluster,
                    "condition": condition,
                }
            )
            .drop_duplicates("cell_identity")
        )
        df = df.merge(meta, on="cell_identity", how="left")

        # Normalise expression for coloring
        df["logMeanExpression"] = np.log2(1.0 + df["meanExpr"])
        df["maxMeanExpression"] = df.groupby("gene")["meanExpr"].transform("max")
        df["relMeanExpression"] = np.where(
            df["maxMeanExpression"] > 0,
            df["meanExpr"] / df["maxMeanExpression"],
            0.0,
        )

        # Order axes nicely
        df["gene"] = pd.Categorical(df["gene"], categories=genes, ordered=True)
        df["cell_identity"] = pd.Categorical(
            df["cell_identity"],
            categories=group_names,
            ordered=True,
        )

        return df.sort_values(["gene", "cell_identity"]).reset_index(drop=True)

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