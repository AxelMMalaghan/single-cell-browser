import anndata as ad
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import pytest
from typing import cast

from sc_browser.core.dataset import Dataset
from sc_browser.core.filter_state import FilterState
from sc_browser.views.volcano_plot_view import VolcanoPlotView


def _make_dataset_for_volcano():
    """
    Small AnnData with:
    - 4 cells
    - clusters: A, A, B, B
    - conditions: ctrl, treated, ctrl, treated
    """
    obs = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3", "c4"],
            "cluster": ["A", "A", "B", "B"],
            "condition": ["ctrl", "treated", "ctrl", "treated"],
            "sample": ["s1", "s1", "s2", "s2"],
        },
        index=pd.Index(["c1", "c2", "c3", "c4"]),
    )

    var = pd.DataFrame(index=pd.Index(["g1", "g2", "g3"]))
    X = np.arange(obs.shape[0] * var.shape[0]).reshape(obs.shape[0], var.shape[0])

    adata = ad.AnnData(X=X, obs=obs, var=var)

    ds = Dataset(
        name="VolcanoDataset",
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


def _make_state_for_volcano(
    conditions=None, clusters=None, samples=None
) -> FilterState:
    """
    Minimal FilterState for VolcanoPlotView.
    """
    conditions = conditions or []
    clusters = clusters or []
    samples = samples or []

    state = FilterState(
        dataset_name="dataset_name",
        view_id="subset",
        clusters=clusters,
        conditions=conditions,
        genes=[],
    )

    # VolcanoPlotView looks at:
    # - state.conditions
    # - state.clusters
    # - state.samples (if present)
    if hasattr(state, "samples"):
        state.samples = samples

    return state


class DummyDEResult:
    def __init__(self, table: pd.DataFrame):
        self.table = table


@pytest.fixture
def fake_de_table():
    """
    Simple DE result with 3 genes:
    - g1: strong up (log2FC=2, p=0.001)
    - g2: not significant (log2FC=0.5, p=0.5)
    - g3: strong down (log2FC=-2, p=0.0001)
    """
    return pd.DataFrame(
        {
            "gene": ["g1", "g2", "g3"],
            "log2FC": [2.0, 0.5, -2.0],
            "pvalue": [0.001, 0.5, 0.0001],
            "adj_pvalue": [0.001, 0.5, 0.0001],
        }
    )


def test_volcano_compute_data_basic(monkeypatch, fake_de_table):
    ds = _make_dataset_for_volcano()
    view = VolcanoPlotView(dataset=ds)

    state = _make_state_for_volcano(
        conditions=["ctrl", "treated"],
        clusters=[],
        samples=[],
    )

    # Patch run_de used inside VolcanoPlotView to return our dummy result
    import sc_browser.views.volcano_plot_view as volcano_mod

    def fake_run_de(config):
        # We don't care about DEConfig internals here; just return the table
        return DummyDEResult(fake_de_table.copy())

    monkeypatch.setattr(volcano_mod, "run_de", fake_run_de)

    df = view.compute_data(state)

    # Basic structure
    assert isinstance(df, pd.DataFrame)
    assert not df.empty

    # Required columns
    for col in ["gene", "log2FC", "neg_log10_pvalue", "significance"]:
        assert col in df.columns

    # Thresholds from implementation
    log_fc_thr = df.attrs["log_fc_threshold"]
    p_thr = df.attrs["pval_threshold"]
    assert log_fc_thr == 1.0
    assert p_thr == 0.05

    # Comparison string
    assert df.attrs["comparison"] in (
        "ctrl vs treated",
        "treated vs ctrl",
        "ctrl vs rest",
        "treated vs rest",
    )

    # Check significance labels
    sig_by_gene = dict(zip(df["gene"], df["significance"], strict=False))

    # g1: log2FC=2, p=0.001 -> Upregulated
    assert sig_by_gene["g1"] == "Upregulated"

    # g3: log2FC=-2, p=0.0001 -> Downregulated
    assert sig_by_gene["g3"] == "Downregulated"

    # g2: log2FC=0.5, p=0.5 -> Not Significant
    assert sig_by_gene["g2"] == "Not Significant"


def test_volcano_compute_data_empty_de_result(monkeypatch):
    ds = _make_dataset_for_volcano()
    view = VolcanoPlotView(dataset=ds)
    state = _make_state_for_volcano(conditions=["ctrl", "treated"])

    import sc_browser.views.volcano_plot_view as volcano_mod

    def fake_run_de(config):
        return DummyDEResult(pd.DataFrame())

    monkeypatch.setattr(volcano_mod, "run_de", fake_run_de)

    df = view.compute_data(state)
    assert isinstance(df, pd.DataFrame)
    assert df.empty


def test_volcano_render_figure_basic(monkeypatch, fake_de_table):
    ds = _make_dataset_for_volcano()
    view = VolcanoPlotView(dataset=ds)
    state = _make_state_for_volcano(conditions=["ctrl", "treated"])

    import sc_browser.views.volcano_plot_view as volcano_mod

    def fake_run_de(config):
        return DummyDEResult(fake_de_table.copy())

    monkeypatch.setattr(volcano_mod, "run_de", fake_run_de)

    df = view.compute_data(state)
    fig = cast(go.Figure, view.render_figure(df, state))

    assert isinstance(fig, go.Figure)
    # axes should be labelled
    assert fig.layout.xaxis.title.text == "log2(Fold Change)"
    assert fig.layout.yaxis.title.text == "-log10(p-value)"
    # at least one trace
    traces = list(fig.data)
    assert len(traces) >= 1


def test_volcano_render_figure_empty():
    ds = _make_dataset_for_volcano()
    view = VolcanoPlotView(dataset=ds)
    state = _make_state_for_volcano(conditions=["ctrl", "treated"])

    empty_df = pd.DataFrame()
    fig = view.render_figure(empty_df, state)

    assert isinstance(fig, go.Figure)
    assert (fig.layout.title.text or "") == "No data to show"
