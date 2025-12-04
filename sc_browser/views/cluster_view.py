from __future__ import annotations

import pandas as pd
import plotly.graph_objs as go

from sc_browser.core.base_view import BaseView
from sc_browser.core.filter_state import FilterState, FilterProfile


class ClusterView(BaseView):
    id = "cluster"
    label = "Clusters"
    filter_profile = FilterProfile(embedding=True)

    # -----------------------------------------------------------
    # compute_data
    # -----------------------------------------------------------
    def compute_data(self, state: FilterState) -> pd.DataFrame:
        ds = self.dataset.subset(
            clusters=state.clusters or None,
            conditions=state.conditions or None,
            samples=state.samples or None,
        )

        adata = ds.adata

        if adata.n_obs == 0:
            return pd.DataFrame()

        # Determine embedding
        emb_key = state.embedding or ds.embedding_key
        if emb_key not in ds.adata.obsm:
            raise ValueError(f"Embedding '{emb_key}' not found in .obsm")

        # Unified Dataset API
        coords = ds.get_embedding_matrix(emb_key)
        labels = ds.get_embedding_labels(emb_key)

        # Build dataframe
        df = pd.DataFrame(
            {
                "x": coords[:, 0],
                "y": coords[:, 1],
            },
            index=adata.obs_names,
        )

        # Optional z dimension
        if coords.shape[1] >= 3:
            df["z"] = coords[:, 2]

        # Store axis labels for renderer
        df.attrs["embedding_labels"] = labels

        # Metadata
        df["cluster"] = ds.clusters.astype(str).values
        df["condition"] = (
            ds.conditions.values if ds.conditions is not None else "all"
        )

        return df

    # -----------------------------------------------------------
    # render_figure
    # -----------------------------------------------------------
    def render_figure(self, data: pd.DataFrame, state: FilterState) -> go.Figure:
        if data.empty:
            fig = go.Figure()
            fig.add_annotation(
                text="No cells after filtering",
                showarrow=False,
                xref="paper",
                yref="paper",
                x=0.5,
                y=0.5,
            )
            return fig

        # 3D only if requested AND z-dim exists
        if state.is_3d and "z" in data.columns:
            return self._render_3d(data)
        else:
            return self._render_2d(data)

    # -----------------------------------------------------------
    # 2D renderer
    # -----------------------------------------------------------
    def _render_2d(self, data: pd.DataFrame) -> go.Figure:
        labels = data.attrs.get("embedding_labels", ["Dim 1", "Dim 2"])

        fig = go.Figure()

        for cluster, dfc in data.groupby("cluster"):
            fig.add_trace(
                go.Scattergl(
                    x=dfc["x"],
                    y=dfc["y"],
                    mode="markers",
                    name=str(cluster),
                    marker={"size": 5, "opacity": 0.7},
                )
            )

        fig.update_layout(
            xaxis_title=labels[0],
            yaxis_title=labels[1],
            legend_title="Cluster",
            margin=dict(l=40, r=40, t=40, b=40),
        )

        return fig

    # -----------------------------------------------------------
    # 3D renderer
    # -----------------------------------------------------------
    def _render_3d(self, data: pd.DataFrame) -> go.Figure:
        labels = data.attrs.get(
            "embedding_labels",
            ["Dim 1", "Dim 2", "Dim 3"]
        )

        fig = go.Figure()

        for cluster, dfc in data.groupby("cluster"):
            fig.add_trace(
                go.Scatter3d(
                    x=dfc["x"],
                    y=dfc["y"],
                    z=dfc["z"],
                    mode="markers",
                    name=str(cluster),
                    marker=dict(size=3, opacity=0.7),
                )
            )

        fig.update_layout(
            scene=dict(
                xaxis_title=labels[0],
                yaxis_title=labels[1],
                zaxis_title=labels[2] if len(labels) > 2 else "Dim 3",
            ),
            legend_title="Cluster",
            margin=dict(l=0, r=0, t=40, b=0),
        )

        return fig