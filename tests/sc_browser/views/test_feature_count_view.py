import anndata as ad
import numpy as np
import pandas as pd
import plotly.graph_objs as go
from typing import cast

from sc_browser.core.dataset import Dataset
from sc_browser.core.filter_state import FilterState
from sc_browser.views.feature_count_view import FeatureCountView


def _make_dataset_for_feature_count():
    """
    Tiny AnnData with:
    - 3 cells
    - 3 genes
    - clusters A, A, B
    - conditions x, y, x
    """
    obs = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3"],
            "cluster": ["A", "A", "B"],
            "condition": ["x", "y", "x"],
            "sample": ["s1", "s1", "s2"],
        },
        index=pd.Index(["c1", "c2", "c3"]),
    )

    var = pd.DataFrame(index=pd.Index(["g1", "g2", "g3"]))
    # rows: cells, cols: genes
    X = np.array(
        [
            [1, 0, 2],  # c1: counts=3,  features>0 = 2
            [0, 5, 0],  # c2: counts=5,  features>0 = 1
            [3, 4, 0],  # c3: counts=7,  features>0 = 2
        ]
    )

    adata = ad.AnnData(X=X, obs=obs, var=var)

    ds = Dataset(
        name="FeatureDataset",
        group="TestGroup",
        adata=adata,
        cluster_key="cluster",
        condition_key="condition",
        embedding_key=None,
        obs_columns={
            "cell_id": "cell_id",
            "cluster": "cluster",
            "condition": "condition",
            "sample": "sample",
        },
        file_path=None,
    )
    return ds


def _make_state(
    split_by_condition: bool = False, clusters=None, conditions=None
) -> FilterState:
    clusters = clusters or []
    conditions = conditions or []

    state = FilterState(
        dataset_name="dataset_name",
        view_id="subset",
        clusters=clusters,
        conditions=conditions,
        genes=[],
    )

    # FeatureCountView expects split_by_condition
    state.split_by_condition = split_by_condition

    return state


def test_feature_count_compute_data_basic():
    ds = _make_dataset_for_feature_count()
    view = FeatureCountView(dataset=ds)
    state = _make_state(split_by_condition=False)

    df = view.compute_data(state)

    assert isinstance(df, pd.DataFrame)
    assert len(df) == ds.adata.n_obs
    assert list(df.index) == list(ds.adata.obs_names)

    # required columns
    for col in ["n_counts", "n_features", "cluster", "condition", "group"]:
        assert col in df.columns

    # counts
    assert list(df["n_counts"]) == [3, 5, 7]
    assert list(df["n_features"]) == [2, 1, 2]

    # clusters / conditions preserved
    assert list(df["cluster"]) == ["A", "A", "B"]
    assert list(df["condition"]) == ["x", "y", "x"]

    # when split_by_condition=False, group == cluster
    assert list(df["group"]) == list(df["cluster"])


def test_feature_count_compute_data_split_by_condition():
    ds = _make_dataset_for_feature_count()
    view = FeatureCountView(dataset=ds)
    state = _make_state(split_by_condition=True)

    df = view.compute_data(state)

    # group should be cluster_condition
    # cells: (A,x), (A,y), (B,x) -> "A_x", "A_y", "B_x"
    assert list(df["group"]) == ["A_x", "A_y", "B_x"]


def test_feature_count_compute_data_with_filtering():
    ds = _make_dataset_for_feature_count()
    view = FeatureCountView(dataset=ds)

    # Filter to cluster B only -> c3
    state = _make_state(clusters=["B"])
    df = view.compute_data(state)

    assert list(df.index) == ["c3"]
    assert list(df["n_counts"]) == [7]
    assert list(df["cluster"]) == ["B"]


def test_feature_count_render_figure_basic():
    ds = _make_dataset_for_feature_count()
    view = FeatureCountView(dataset=ds)
    state = _make_state()

    df = view.compute_data(state)
    fig = cast(go.Figure, view.render_figure(df, state))

    assert isinstance(fig, go.Figure)

    # expect at least one trace
    traces = list(fig.data)
    assert len(traces) >= 1

    # axis titles should be set
    assert fig.layout.xaxis.title.text == "Total counts per cell"
    assert fig.layout.yaxis.title.text == "Number of features per cell"


def test_feature_count_render_figure_empty():
    ds = _make_dataset_for_feature_count()
    view = FeatureCountView(dataset=ds)

    # Force empty by picking a non-existent cluster
    state = _make_state(clusters=["Z"])
    df = view.compute_data(state)

    assert df.empty

    fig = view.render_figure(df, state)
    assert isinstance(fig, go.Figure)
