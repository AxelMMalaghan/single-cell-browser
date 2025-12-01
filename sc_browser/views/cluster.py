from __future__ import annotations

import pandas as pd
import plotly.express as px
from plotly.graph_objs import Figure

from ..core.state import FilterState
from ..core.base_view import BaseView


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
        ds = self.dataset.subset(
            clusters = state.clusters or None,
            conditions = state.conditions or None,
        )

        embedding = ds.adata.obsm[ds.embedding_key]

        df = pd.DataFrame(embedding[:, :2], columns=["dim1", "dim2"])
        df["cluster"] = ds.clusters.to_numpy()
        df["condition"] = ds.conditions.to_numpy()
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