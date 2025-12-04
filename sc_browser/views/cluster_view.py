from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
from plotly.graph_objs import Figure

from sc_browser.core.dataset import Dataset
from sc_browser.core.state import FilterState, FilterProfile
from sc_browser.core.base_view import BaseView


class ClusterView(BaseView):
    """
    UMAP/TSNE/PCA /.... cluster plot

    - X/Y from embedding
    - Color by cluster
    - Optional facet by condition
    """

    id = "cluster"
    label = "Cluster Plot"

    # This view cares about clusters, conditions, samples, cell types, and embedding;
    # no gene widget needed.
    filter_profile = FilterProfile(
        clusters=True,
        conditions=True,
        samples=True,
        cell_types=True,
        genes=False,
        embedding=True,
    )

    def compute_data(self, state: FilterState) -> pd.DataFrame:
        ds: Dataset = self.dataset

        # Apply filters using Dataset abstraction (this hits the subset cache)
        ds_sub = ds.subset(
            clusters=state.clusters or None,
            conditions=state.conditions or None,
            samples=state.samples or None,
            cell_types=state.cell_types or None,
        )

        # No cells after filtering → return empty frame
        if ds_sub.adata.n_obs == 0:
            return pd.DataFrame(columns=["dim1", "dim2", "cluster", "condition"])

        # Decide which embedding key to use:
        #  - use the widget selection if set
        #  - otherwise fall back to the dataset's configured embedding_key
        emb_key = state.embedding or ds_sub.embedding_key
        if emb_key not in ds_sub.adata.obsm:
            raise ValueError(
                f"Selected embedding '{emb_key}' not found in adata.obsm. "
                f"Available keys: {list(ds_sub.adata.obsm.keys())}"
            )

        emb = ds_sub.adata.obsm[emb_key]

        # Handle both DataFrame and array-like obsm storage
        if isinstance(emb, pd.DataFrame):
            arr = emb.to_numpy()
        else:
            arr = np.asarray(emb)

        if arr.ndim != 2 or arr.shape[1] < 2:
            raise ValueError(
                f"Embedding '{emb_key}' must be a 2D array with at least 2 columns, "
                f"got shape {arr.shape}"
            )

        # Build embedding DataFrame with dim1/dim2
        emb_df = pd.DataFrame(
            arr[:, :2],
            index=ds_sub.adata.obs_names,
            columns=["dim1", "dim2"],
        )

        # Attach cluster labels (pre-normalised in Dataset)
        emb_df["cluster"] = ds_sub.clusters.astype(str).values

        # Attach condition only if available
        if ds_sub.condition_key is not None and ds_sub.conditions is not None:
            emb_df["condition"] = ds_sub.conditions.astype(str).values
        else:
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

        facet_col = (
            "condition"
            if state.split_by_condition and "condition" in data.columns
            else None
        )

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