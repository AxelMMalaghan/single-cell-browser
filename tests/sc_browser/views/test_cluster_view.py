import numpy as np
import pandas as pd
import anndata as ad
import plotly.graph_objs as go

from sc_browser.core.dataset import Dataset
from sc_browser.core.filter_state import FilterState
from sc_browser.views.cluster_view import ClusterView


def _make_dataset_with_embedding():
    """
    Tiny AnnData with:
    - 3 cells
    - 2 clusters (A, B)
    - 2 conditions (x, y)
    - 2D embedding in .obsm["X_umap"]
    """
    obs = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3"],
            "cluster": ["A", "A", "B"],
            "condition": ["x", "y", "x"],
        },
        index=["c1", "c2", "c3"],
    )

    var = pd.DataFrame(index=["g1", "g2"])
    X = np.array(
        [
            [1.0, 2.0],
            [3.0, 4.0],
            [5.0, 6.0],
        ]
    )

    adata = ad.AnnData(X=X, obs=obs, var=var)

    # 2D embedding: simple coordinates
    emb = np.array(
        [
            [0.0, 0.0],  # c1
            [1.0, 0.0],  # c2
            [0.0, 1.0],  # c3
        ]
    )
    adata.obsm["X_umap"] = emb

    ds = Dataset(
        name="TestDataset",
        group="TestGroup",
        adata=adata,
        cluster_key="cluster",
        condition_key="condition",
        embedding_key="X_umap",
        obs_columns={
            "cell_id": "cell_id",
            "cluster": "cluster",
            "condition": "condition",
        },
        file_path=None,
    )

    return ds


def _make_default_state() -> FilterState:
    """
    Minimal FilterState for ClusterView tests.
    We only care about embedding + 2D mode; other fields can use defaults.
    """
    state = FilterState(
        dataset_name="dataset_name",
        view_id="subset",
        clusters=[],
        conditions=[],
        genes=[],
    )

    # Attributes ClusterView expects but FilterState may not set in __init__
    state.embedding = None  # use dataset.default embedding_key ("X_umap")
    state.dim_x = None
    state.dim_y = None
    state.dim_z = None
    state.is_3d = False

    return state


def test_cluster_view_compute_data_basic():
    ds = _make_dataset_with_embedding()
    view = ClusterView(dataset=ds)  # BaseView should accept dataset=...

    state = _make_default_state()

    data = view.compute_data(state)

    # DataFrame shape and index
    assert isinstance(data, pd.DataFrame)
    assert len(data) == ds.adata.n_obs
    assert list(data.index) == list(ds.adata.obs_names)

    # Required columns
    for col in ["x", "y", "cluster", "condition"]:
        assert col in data.columns

    # Clusters/conditions preserved
    assert list(data["cluster"]) == ["A", "A", "B"]
    assert list(data["condition"]) == ["x", "y", "x"]

    # Axis labels derived from "X_umap"
    # X_umap -> "UMAP 1", "UMAP 2"
    assert data.attrs["x_label"] == "UMAP 1"
    assert data.attrs["y_label"] == "UMAP 2"


def test_cluster_view_render_figure_2d():
    ds = _make_dataset_with_embedding()
    view = ClusterView(dataset=ds)
    state = _make_default_state()

    data = view.compute_data(state)
    fig = view.render_figure(data, state)

    assert isinstance(fig, go.Figure)

    # One trace per cluster (A, B) -> 2 traces
    clusters = sorted(data["cluster"].unique())
    assert len(fig.data) == len(clusters)

    # Expect scattergl traces for 2D mode
    for trace in fig.data:
        assert trace.type.lower() in ("scattergl", "scatter")

    # Axis titles from attrs
    layout = fig.layout
    assert layout.xaxis.title.text == data.attrs["x_label"]
    assert layout.yaxis.title.text == data.attrs["y_label"]


def test_cluster_view_render_figure_empty():
    ds = _make_dataset_with_embedding()
    view = ClusterView(dataset=ds)
    state = _make_default_state()

    # Simulate no data after filtering
    empty_df = pd.DataFrame()
    fig = view.render_figure(empty_df, state)

    # Should still return a valid Figure (BaseView.empty_figure)
    assert isinstance(fig, go.Figure)
