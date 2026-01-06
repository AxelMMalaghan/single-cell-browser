import anndata as ad
import numpy as np
import pandas as pd
import plotly.graph_objs as go
from typing import cast

from sc_browser.core.dataset import Dataset
from sc_browser.core.filter_state import FilterState
from sc_browser.views.heat_map_view import HeatmapView


def _make_dataset_for_heatmap():
    """
    Tiny AnnData with:
    - 4 cells
    - 3 genes
    - clusters A, A, B, B
    - conditions x, y, x, y
    - samples s1,s1,s2,s2
    """
    obs = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3", "c4"],
            "cluster": ["A", "A", "B", "B"],
            "condition": ["x", "y", "x", "y"],
            "sample": ["s1", "s1", "s2", "s2"],
        },
        index=pd.Index(["c1", "c2", "c3", "c4"]),
    )

    var = pd.DataFrame(index=pd.Index(["g1", "g2", "g3"]))

    # rows: cells, cols: genes
    # g1: [1, 3, 5, 7]
    # g2: [2, 4, 6, 8]
    # g3: [0, 1, 0, 1]
    X = np.array(
        [
            [1, 2, 0],
            [3, 4, 1],
            [5, 6, 0],
            [7, 8, 1],
        ]
    )

    adata = ad.AnnData(X=X, obs=obs, var=var)

    ds = Dataset(
        name="HeatmapDataset",
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
    *,
    genes=None,
    split_by_condition: bool = False,
    clusters=None,
    conditions=None,
    samples=None,
    cell_types=None,
    color_scale: str = "viridis",
) -> FilterState:
    genes = genes or []
    clusters = clusters or []
    conditions = conditions or []
    samples = samples or []
    cell_types = cell_types or []

    state = FilterState(
        dataset_name="HeatmapDataset",
        view_id="heatmap",
        genes=genes,
        clusters=clusters,
        conditions=conditions,
        samples=samples,
        cell_types=cell_types,
        color_scale=color_scale,
    )

    state.split_by_condition = split_by_condition
    return state


def test_heatmap_compute_data_no_genes_returns_empty_df():
    ds = _make_dataset_for_heatmap()
    view = HeatmapView(dataset=ds)
    state = _make_state(genes=[])

    df = view.compute_data(state)

    assert isinstance(df, pd.DataFrame)
    assert df.empty


def test_heatmap_compute_data_filters_genes_to_present_only():
    ds = _make_dataset_for_heatmap()
    view = HeatmapView(dataset=ds)

    # g999 not present, should be dropped, g1 should survive
    state = _make_state(genes=["g1", "g999"], split_by_condition=False)
    df = view.compute_data(state)

    assert not df.empty
    assert set(df["gene"].unique()) == {"g1"}


def test_heatmap_compute_data_split_by_condition_groups_correctly():
    ds = _make_dataset_for_heatmap()
    view = HeatmapView(dataset=ds)

    # With split_by_condition=True, groups should be:
    # A_x (c1), A_y (c2), B_x (c3), B_y (c4)
    state = _make_state(genes=["g1", "g2"], split_by_condition=True)
    df = view.compute_data(state)

    assert not df.empty

    # Required columns
    for col in ["group", "gene", "mean_expression", "log_mean_expression"]:
        assert col in df.columns

    # Groups produced
    assert set(df["group"].unique()) == {"A_x", "A_y", "B_x", "B_y"}

    # One row per (group, gene) in this tiny dataset (each group has 1 cell)
    # => 4 groups * 2 genes = 8 rows
    assert len(df) == 8

    # Check one concrete value: A_x, g1 mean = 1 -> log1p(1)
    row = cast(pd.DataFrame, df[(df["group"] == "A_x") & (df["gene"] == "g1")])
    assert len(row) == 1
    assert float(row["mean_expression"].iloc[0]) == 1.0
    assert np.isclose(float(row["log_mean_expression"].iloc[0]), np.log1p(1.0))


def test_heatmap_compute_data_no_split_groups_by_cluster_only():
    ds = _make_dataset_for_heatmap()
    view = HeatmapView(dataset=ds)

    state = _make_state(genes=["g1", "g2"], split_by_condition=False)
    df = view.compute_data(state)

    assert not df.empty
    assert set(df["group"].unique()) == {"A", "B"}

    # Now groups are A and B (each has 2 cells), genes 2 => 4 rows
    assert len(df) == 4

    # Check mean for group A, gene g1: cells c1=1, c2=3 => mean 2
    row = cast(pd.DataFrame, df[(df["group"] == "A") & (df["gene"] == "g1")])
    assert len(row) == 1
    assert float(row["mean_expression"].iloc[0]) == 2.0
    assert np.isclose(float(row["log_mean_expression"].iloc[0]), np.log1p(2.0))


def test_heatmap_compute_data_with_filtering_clusters():
    ds = _make_dataset_for_heatmap()
    view = HeatmapView(dataset=ds)

    # Filter to cluster B only -> cells c3,c4
    state = _make_state(genes=["g1"], clusters=["B"], split_by_condition=False)
    df = view.compute_data(state)

    assert not df.empty
    assert set(df["group"].unique()) == {"B"}

    # Mean for B, g1: c3=5, c4=7 => mean 6
    row = cast(pd.DataFrame, df[(df["group"] == "B") & (df["gene"] == "g1")])
    assert float(row["mean_expression"].iloc[0]) == 6.0
    assert np.isclose(float(row["log_mean_expression"].iloc[0]), np.log1p(6.0))


def test_heatmap_render_figure_basic():
    ds = _make_dataset_for_heatmap()
    view = HeatmapView(dataset=ds)
    state = _make_state(
        genes=["g1", "g2"], split_by_condition=False, color_scale="viridis"
    )

    df = view.compute_data(state)
    fig = cast(go.Figure, view.render_figure(df, state))

    assert isinstance(fig, go.Figure)
    traces = list(fig.data)
    assert len(traces) >= 1

    # px.imshow sets a coloraxis title via colorbar title; check it exists
    # (plotly can store it in different places depending on version)
    assert fig.layout.height == 600
    assert fig.layout.margin.l == 40
    assert fig.layout.margin.r == 40
    assert fig.layout.margin.t == 40
    assert fig.layout.margin.b == 40


def test_heatmap_render_figure_empty():
    ds = _make_dataset_for_heatmap()
    view = HeatmapView(dataset=ds)

    # Force empty: no genes
    state = _make_state(genes=[])
    df = view.compute_data(state)
    assert df.empty

    fig = view.render_figure(df, state)
    assert isinstance(fig, go.Figure)
