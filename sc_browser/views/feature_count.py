from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.graph_objs import Figure

from sc_browser.core.dataset import Dataset
from sc_browser.core.view_base import BaseView
from sc_browser.core.state import FilterState

class FeatureCountView(BaseView):
    """
    QC plot: per cell total counts vs number of detected features
    """

    id = "feature_count"
    label = "Feature v. Count"

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

        if adata_sub.n_obs == 0:
            return pd.DataFrame()

        # --- count matrix (cells x genes) ---
        X = adata_sub.X
        if hasattr(X, "toarray"):
            X = X.toarray()

        # per-cell total counts and detected features
        n_counts = X.sum(axis=1)
        n_features = (X > 0).sum(axis=1)

        cluster = obs_sub[self.dataset.cluster_key].astype(str)
        condition = obs_sub[self.dataset.condition_key].astype(str)

        if state.split_by_condition:
            group = cluster + "-" + condition
        else:
            group = cluster

        df = pd.DataFrame(
            {
                "n_counts": np.asarray(n_counts).ravel(),
                "n_features": np.asarray(n_features).ravel(),
                "cluster": cluster.values,
                "condition": condition.values,
                "group": group.values,
            },
            index=obs_sub.index,
        )
        return df

    def render_figure(self, data: pd.DataFrame, state: FilterState) -> Any:

        if data is None or data.empty:
            fig = go.Figure()
            fig.update_layout(
                title="No cells after filtering - adjust cluster/condition filters",
                xaxis={"visible": False},
                yaxis={"visible": False},
            )
            return fig

        fig = px.scatter(
            data,
            x="n_counts",
            y="n_features",
            color="group",
            hover_data={"cluster": True, "condition": True, "n_counts": True, "n_features": True},
        )

        fig.update_layout(
            height=600,
            margin=dict(l=40, r=40, t=40, b=40),
            xaxis_title="Total counts per cell",
            yaxis_title="Number of features per cell",
            legend_title="Group"
        )
        return fig

