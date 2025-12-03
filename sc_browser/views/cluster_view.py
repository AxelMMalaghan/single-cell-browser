from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
from plotly.graph_objs import Figure

from sc_browser.core import Dataset
from sc_browser.core.state import FilterState
from sc_browser.core.base_view import BaseView


class ClusterView(BaseView):
    """
    UMAP/TSNE cluster plot

    - X/Y from embedding
    - Color by cluster
    - Optional facet by condition
    """

    id = "cluster"
    label = "Cluster Plot"

    def compute_data(self, state: FilterState) -> pd.DataFrame:
        ds: Dataset = self.dataset

        # Apply filters using Dataset abstraction
        ds_sub = ds.subset(
            clusters=state.clusters or None,
            conditions=state.conditions or None,
        )
        adata_sub = ds_sub.adata

        # Get embedding as array (2D)
        emb = ds_sub.adata.obsm[ds_sub.embedding_key]
        if isinstance(emb, pd.DataFrame):
            emb_arr = emb.values
        else:
            emb_arr = np.asarray(emb)

        df = pd.DataFrame(
            emb_arr[:, :2],
            columns=["dim1", "dim2"],
            index=adata_sub.obs_names,
        )

        # 🔹 Add cluster / condition columns for Plotly
        df["cluster"] = adata_sub.obs[ds_sub.cluster_key].astype(str).values

        # Optional: only if you use it elsewhere
        if ds_sub.condition_key:
            df["condition"] = adata_sub.obs[ds_sub.condition_key].astype(str).values

        return df



    def render_figure(self, data: pd.DataFrame, state: FilterState) -> Figure:

        facet_col = "condition" if state.split_by_condition else None

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