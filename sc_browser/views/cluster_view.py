from __future__ import annotations

import pandas as pd
import plotly.graph_objs as go

from sc_browser.core.base_view import BaseView
from sc_browser.core.filter_state import FilterState
from sc_browser.core.filter_profile import FilterProfile


class ClusterView(BaseView):
    id = "cluster"
    label = "Clusters"
    filter_profile = FilterProfile(
        clusters=True,
        conditions=True,
        samples=True,
        cell_types=True,
        genes=False,
        embedding=True,
        split_by_condition=True,
        is_3d=True,
    )

    # -----------------------------------------------------------
    # Derive human-readable axis names from embedding key
    # -----------------------------------------------------------
    def _derive_axis_labels(self, emb_key: str, n_dims: int) -> list[str]:
        name = emb_key
        if name.lower().startswith("x_"):
            name = name[2:]  # strip leading X_
        pretty = name.replace("_", " ").upper()
        return [f"{pretty} {i+1}" for i in range(n_dims)]

    # -----------------------------------------------------------
    # Compute dataframe for rendering
    # -----------------------------------------------------------
    def compute_data(self, state: FilterState) -> pd.DataFrame:
        # Centralised filtering: let Dataset handle cluster/condition/sample/cell-type
        ds = self.filtered_dataset(state)
        adata = ds.adata

        if adata.n_obs == 0:
            return pd.DataFrame()

        # Choose embedding: per-state override or dataset default
        emb_key = state.embedding or ds.embedding_key
        if emb_key not in ds.adata.obsm:
            raise ValueError(f"Embedding '{emb_key}' not found in .obsm")

        coords = ds.get_embedding_matrix(emb_key)
        labels = self._derive_axis_labels(emb_key, coords.shape[1])

        df = pd.DataFrame(index=adata.obs_names)

        # Safe dimension defaults
        dim_x = state.dim_x if state.dim_x is not None else 0
        dim_y = state.dim_y if state.dim_y is not None else 1
        dim_z = state.dim_z if state.dim_z is not None else 2

        df["x"] = coords[:, dim_x]
        df["y"] = coords[:, dim_y]

        if state.is_3d and coords.shape[1] > dim_z:
            df["z"] = coords[:, dim_z]

        # Store embedding labels for axis titles
        df.attrs["embedding_labels"] = labels
        df.attrs["x_label"] = labels[dim_x]
        df.attrs["y_label"] = labels[dim_y]
        df.attrs["z_label"] = labels[dim_z] if len(labels) > dim_z else None

        # Metadata
        if ds.clusters is not None:
            df["cluster"] = ds.clusters.astype(str).values
        else:
            df["cluster"] = "NA"

        if ds.conditions is not None:
            df["condition"] = ds.conditions.values
        else:
            df["condition"] = "all"

        return df

    # -----------------------------------------------------------
    # Render entrypoint
    # -----------------------------------------------------------
    def render_figure(self, data: pd.DataFrame, state: FilterState) -> go.Figure:

        if data is None or data.empty:
            return self.empty_figure("No data to show")

        if state.is_3d and "z" in data.columns:
            return self._render_3d(data)
        return self._render_2d(data)

    # -----------------------------------------------------------
    # 2D Renderer
    # -----------------------------------------------------------
    def _render_2d(self, data: pd.DataFrame) -> go.Figure:
        x_label = data.attrs["x_label"]
        y_label = data.attrs["y_label"]

        fig = go.Figure()
        for cluster, dfc in data.groupby("cluster"):
            fig.add_trace(
                go.Scattergl(
                    x=dfc["x"],
                    y=dfc["y"],
                    mode="markers",
                    name=str(cluster),
                    marker=dict(size=5, opacity=0.7),
                )
            )

        fig.update_layout(
            xaxis_title=x_label,
            yaxis_title=y_label,
            legend_title="Cluster",
            margin=dict(l=40, r=40, t=40, b=40),
        )
        return fig

    # -----------------------------------------------------------
    # 3D Renderer
    # -----------------------------------------------------------
    def _render_3d(self, data: pd.DataFrame) -> go.Figure:
        x_label = data.attrs["x_label"]
        y_label = data.attrs["y_label"]
        z_label = data.attrs["z_label"]

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
                xaxis_title=x_label,
                yaxis_title=y_label,
                zaxis_title=z_label,
            ),
            legend_title="Cluster",
            margin=dict(l=0, r=0, t=40, b=0),
        )
        return fig
