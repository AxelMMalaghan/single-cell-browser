from typing import cast

import anndata as ad
import numpy as np
import pandas as pd
import plotly.graph_objs as go

from sc_browser.core.dataset import Dataset
from sc_browser.core.filter_state import FilterState
from sc_browser.views.dataset_summary_view import DatasetSummary


def _make_dataset_for_summary():
    """
    Tiny AnnData with:
    - 4 cells
    - 2 genes
    - clusters: A, A, B, B
    - conditions: x, y, x, y
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
    X = np.arange(obs.shape[0] * var.shape[0]).reshape(obs.shape[0], var.shape[0])

    adata = ad.AnnData(X=X, obs=obs, var=var)

    ds = Dataset(
        name="SummaryDataset",
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


def _make_state() -> FilterState:
    """
    Minimal FilterState for DatasetSummary.
    """
    state = FilterState(
        dataset_name="dataset_name",
        view_id="subset",
        clusters=[],  # no filtering
        conditions=[],
        genes=[],
    )
    return state


def test_dataset_summary_compute_data_basic():
    ds = _make_dataset_for_summary()
    view = DatasetSummary(dataset=ds)
    state = _make_state()

    data = view.compute_data(state)

    # compute_data returns a dict, not a DataFrame
    assert isinstance(data, dict)

    # Required keys
    for key in [
        "n_cells",
        "n_genes",
        "obs_schema",
        "cluster_counts",
        "condition_counts",
    ]:
        assert key in data

    assert data["n_cells"] == 4
    assert data["n_genes"] == 2

    obs_schema = data["obs_schema"]
    cluster_counts = data["cluster_counts"]
    condition_counts = data["condition_counts"]

    # obs_schema is a DataFrame with at least 'column', 'dtype', 'n_unique'
    assert isinstance(obs_schema, pd.DataFrame)
    assert {"column", "dtype", "n_unique"}.issubset(obs_schema.columns)

    # Cluster counts: A=2, B=2 (order may vary)
    assert set(cluster_counts["cluster"]) == {"A", "B"}
    counts_by_cluster = dict(
        zip(cluster_counts["cluster"], cluster_counts["count"], strict=False)
    )
    assert counts_by_cluster["A"] == 2
    assert counts_by_cluster["B"] == 2

    # Condition counts: x=2, y=2
    assert set(condition_counts["condition"]) == {"x", "y"}
    counts_by_cond = dict(
        zip(condition_counts["condition"], condition_counts["count"], strict=False)
    )
    assert counts_by_cond["x"] == 2
    assert counts_by_cond["y"] == 2


def test_dataset_summary_render_figure_basic():
    ds = _make_dataset_for_summary()
    view = DatasetSummary(dataset=ds)
    state = _make_state()

    data = view.compute_data(state)
    fig = cast(go.Figure, view.render_figure(data, state))

    assert isinstance(fig, go.Figure)

    # Title should include n_cells and n_genes
    title_text = fig.layout.title.text
    assert "4 cells" in title_text
    assert "2 genes" in title_text

    # We expect at least one bar trace (clusters) and possibly another (conditions)
    traces = list(fig.data)
    assert len(traces) >= 1


def test_dataset_summary_render_figure_empty():
    ds = _make_dataset_for_summary()
    view = DatasetSummary(dataset=ds)
    state = _make_state()

    empty_data = {}

    fig = view.render_figure(empty_data, state)

    assert isinstance(fig, go.Figure)
    assert (fig.layout.title.text or "") == "No data to show"
