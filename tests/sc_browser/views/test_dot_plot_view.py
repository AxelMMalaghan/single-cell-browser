from typing import cast

import anndata as ad
import numpy as np
import pandas as pd
import plotly.graph_objs as go

from sc_browser.core.dataset import Dataset
from sc_browser.core.filter_state import FilterState
from sc_browser.views.dot_plot_view import DotplotView


def _make_dataset_for_dotplot():
    """
    AnnData with:
    - 4 cells, 2 genes
    - clusters: A, A, B, B
    - conditions: x, y, x, y

    Expression:
        c1: g1=1, g2=0
        c2: g1=0, g2=2
        c3: g1=3, g2=0
        c4: g1=0, g2=4
    """

    obs = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3", "c4"],
            "cluster": ["A", "A", "B", "B"],
            "condition": ["x", "y", "x", "y"],
        },
        index=pd.Index(["c1", "c2", "c3", "c4"]),
    )

    var = pd.DataFrame(index=pd.Index(["g1", "g2"]))
    X = np.array(
        [
            [1, 0],  # c1
            [0, 2],  # c2
            [3, 0],  # c3
            [0, 4],  # c4
        ]
    )

    adata = ad.AnnData(X=X, obs=obs, var=var)

    ds = Dataset(
        name="DotplotDataset",
        group="TestGroup",
        adata=adata,
        cluster_key="cluster",
        condition_key="condition",
        embedding_key=None,
        obs_columns={
            "cell_id": "cell_id",
            "cluster": "cluster",
            "condition": "condition",
        },
        file_path=None,
    )
    return ds


def _make_state_for_dotplot(genes=None, split_by_condition=False) -> FilterState:
    genes = genes or []

    state = FilterState(
        dataset_name="dataset_name",
        view_id="subset",
        clusters=[],
        conditions=[],
        genes=genes,
    )

    state.split_by_condition = split_by_condition
    state.color_scale = "viridis"

    return state


def test_dotplot_no_genes_returns_empty_df():
    ds = _make_dataset_for_dotplot()
    view = DotplotView(dataset=ds)

    state = _make_state_for_dotplot(genes=[])

    df = view.compute_data(state)

    assert isinstance(df, pd.DataFrame)
    assert df.empty


def test_dotplot_compute_data_basic():
    ds = _make_dataset_for_dotplot()
    view = DotplotView(dataset=ds)

    state = _make_state_for_dotplot(genes=["g1", "g2"])

    df = view.compute_data(state)

    assert isinstance(df, pd.DataFrame)
    assert not df.empty

    for col in [
        "gene",
        "cell_identity",
        "pctExpressed",
        "relMeanExpression",
        "meanExpr",
        "logMeanExpression",
        "cluster",
        "condition",
    ]:
        assert col in df.columns

    # We expect rows for each (cluster, condition, gene)
    assert set(df["cluster"]) == {"A", "B"}
    assert set(df["condition"]) == {"x", "y"}
    assert set(df["gene"]) == {"g1", "g2"}

    # Check concrete stats for g1
    g1 = df[df["gene"] == "g1"].set_index(["cluster", "condition"])

    # A,x: value 1 -> 100%
    assert g1.loc[("A", "x"), "pctExpressed"] == 100.0
    assert g1.loc[("A", "x"), "meanExpr"] == 1.0

    # A,y: value 0 -> 0%
    assert g1.loc[("A", "y"), "pctExpressed"] == 0.0
    assert g1.loc[("A", "y"), "meanExpr"] == 0.0

    # B,x: value 3 -> 100%
    assert g1.loc[("B", "x"), "pctExpressed"] == 100.0
    assert g1.loc[("B", "x"), "meanExpr"] == 3.0

    # B,y: value 0 -> 0%
    assert g1.loc[("B", "y"), "pctExpressed"] == 0.0
    assert g1.loc[("B", "y"), "meanExpr"] == 0.0


def test_dotplot_split_by_condition_changes_cell_identity():
    ds = _make_dataset_for_dotplot()
    view = DotplotView(dataset=ds)

    state = _make_state_for_dotplot(genes=["g1"], split_by_condition=True)

    df = view.compute_data(state)

    # cell_identity should be "cluster | condition"
    # combos: (A,x), (A,y), (B,x), (B,y)
    expected_identities = {"A | x", "A | y", "B | x", "B | y"}
    assert set(df["cell_identity"].unique()) == expected_identities


def test_dotplot_render_figure_basic():
    ds = _make_dataset_for_dotplot()
    view = DotplotView(dataset=ds)

    state = _make_state_for_dotplot(genes=["g1", "g2"])

    df = view.compute_data(state)
    fig = cast(go.Figure, view.render_figure(df, state))

    assert isinstance(fig, go.Figure)
    # at least one trace in the figure
    traces = list(fig.data)
    assert len(traces) >= 1

    # axes titles should match implementation
    assert fig.layout.xaxis.title.text == "Gene"
    assert fig.layout.yaxis.title.text == "Cluster / Group"


def test_dotplot_render_figure_empty():
    ds = _make_dataset_for_dotplot()
    view = DotplotView(dataset=ds)
    state = _make_state_for_dotplot(genes=[])

    df = view.compute_data(state)
    assert df.empty

    fig = view.render_figure(df, state)

    assert isinstance(fig, go.Figure)
    assert (fig.layout.title.text or "") == "No data to show"
