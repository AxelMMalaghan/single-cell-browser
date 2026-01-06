import anndata as ad
import numpy as np
import pandas as pd
import plotly.graph_objs as go
from typing import cast

from sc_browser.core.dataset import Dataset
from sc_browser.core.filter_state import FilterState
from sc_browser.views.gene_expression_view import ExpressionView


def _make_dataset_for_expression():
    """
    Tiny AnnData with:
    - 3 cells
    - 1 gene 'G1'
    - clusters A, A, B
    - conditions x, y, x
    - 2D embedding in .obsm["X_umap"]
    """
    obs = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3"],
            "cluster": ["A", "A", "B"],
            "condition": ["x", "y", "x"],
        },
        index=pd.Index(["c1", "c2", "c3"]),
    )

    var = pd.DataFrame(index=pd.Index(["G1"]))
    X = np.array([[1.0], [2.0], [0.0]])  # simple expression for G1

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
        name="ExprDataset",
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


def _make_state_for_expression(genes=None, split_by_condition=False) -> FilterState:
    """
    Minimal FilterState for ExpressionView tests.
    We care about: genes, split_by_condition, embedding.
    """
    genes = genes or []

    state = FilterState(
        dataset_name="dataset_name",
        view_id="subset",
        clusters=[],
        conditions=[],
        samples=[] if hasattr(FilterState, "samples") else [],
        genes=genes,
    )

    # Attributes ExpressionView expects but may not be in __init__
    state.embedding = None  # use dataset.embedding_key ("X_umap")
    state.split_by_condition = split_by_condition
    state.color_scale = "viridis"

    return state


def test_expression_view_no_genes_returns_empty_df():
    ds = _make_dataset_for_expression()
    view = ExpressionView(dataset=ds)

    state = _make_state_for_expression(genes=[])

    df = view.compute_data(state)

    assert isinstance(df, pd.DataFrame)
    assert df.empty


def test_expression_view_compute_data_basic():
    ds = _make_dataset_for_expression()
    view = ExpressionView(dataset=ds)

    state = _make_state_for_expression(genes=["G1"])

    df = view.compute_data(state)

    # We should get one row per cell
    assert isinstance(df, pd.DataFrame)
    assert len(df) == ds.adata.n_obs
    assert list(df.index) == list(ds.adata.obs_names)

    # Required columns
    for col in [
        "x",
        "y",
        "expression",
        "log_expression",
        "cluster",
        "condition",
        "group",
        "gene",
    ]:
        assert col in df.columns

    # Gene column should be constant 'G1'
    assert set(df["gene"].unique()) == {"G1"}

    # Group should default to cluster when split_by_condition=False
    assert list(df["group"]) == list(df["cluster"])

    # Clusters / conditions preserved
    assert list(df["cluster"]) == ["A", "A", "B"]
    assert list(df["condition"]) == ["x", "y", "x"]


def test_expression_view_split_by_condition_changes_group():
    ds = _make_dataset_for_expression()
    view = ExpressionView(dataset=ds)

    state = _make_state_for_expression(genes=["G1"], split_by_condition=True)

    df = view.compute_data(state)

    # group should be cluster_condition
    # cells: (A,x), (A,y), (B,x) -> "A_x", "A_y", "B_x"
    assert list(df["group"]) == ["A_x", "A_y", "B_x"]


def test_expression_view_render_figure_basic():
    ds = _make_dataset_for_expression()
    view = ExpressionView(dataset=ds)

    state = _make_state_for_expression(genes=["G1"])

    df = view.compute_data(state)
    fig = cast(go.Figure, view.render_figure(df, state))

    assert isinstance(fig, go.Figure)

    # Expect a scatter-style figure with x/y axes and a color axis
    assert fig.layout.xaxis.title.text is not None
    assert fig.layout.yaxis.title.text is not None
    # Some data points should be present
    traces = list(fig.data)
    assert len(traces) >= 1
