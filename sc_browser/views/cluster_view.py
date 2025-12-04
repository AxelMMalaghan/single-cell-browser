from __future__ import annotations

import pandas as pd
import plotly.express as px
from plotly.graph_objs import Figure

from sc_browser.core.dataset import Dataset
from sc_browser.core.state import FilterState, FilterProfile
from sc_browser.core.base_view import BaseView


class ClusterView(BaseView):
    """
    UMAP/TSNE/PCA cluster plot

    - X/Y from embedding
    - Color by cluster
    - Optional facet by condition
    """

    id = "cluster"
    label = "Cluster Plot"
    filter_profile = FilterProfile()

    def compute_data(self, state: FilterState) -> pd.DataFrame:
        ds: Dataset = self.dataset

        # Apply filters using Dataset abstraction (this hits the subset cache)
        ds_sub = ds.subset(
            clusters=state.clusters or None,
            conditions=state.conditions or None,
            samples=state.samples or None,
            cell_types=getattr(state, "cell_types", None) or None,
        )


        # No cells after filtering → return empty frame
        if ds_sub.adata.n_obs == 0:
            return pd.DataFrame(columns=["dim1", "dim2", "cluster", "condition"])

        # Use Dataset.embedding property – already returns a DataFrame with dim1/dim2
        emb_df = ds_sub.embedding.copy()

        # Attach cluster labels (pre-normalised in Dataset)
        emb_df["cluster"] = ds_sub.clusters.values

        # Attach condition only if available; we only need the column if split_by_condition is used,
        # but it’s cheap enough to add once here.
        if ds_sub.condition_key is not None and ds_sub.conditions is not None:
            emb_df["condition"] = ds_sub.conditions.values
        else:
            # keep column for consistent schema if you prefer
            emb_df["condition"] = None

        return emb_df

    def render_figure(self, data: pd.DataFrame, state: FilterState) -> Figure:
        # If there is no data, return an empty figure with a friendly title
        if data.empty:
            fig = px.scatter(
                title=f"{self.dataset.name} Cluster Plot (no cells after filtering)"
            )
            fig.update_layout(margin=dict(l=40, r=40, t=40, b=40))
            return fig

        facet_col = "condition" if state.split_by_condition and "condition" in data.columns else None

        fig = px.scatter(
            data,
            x="dim1",
            y="dim2",
            color="cluster",
            facet_col=facet_col,
            title=f"{self.dataset.name} Cluster Plot",
        )
        fig.update_traces(marker=dict(size=4))
        fig.update_layout(margin=dict(l=40, r=40, t=40, b=40))
        return fig
