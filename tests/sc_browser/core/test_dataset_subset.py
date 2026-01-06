import numpy as np
import pandas as pd
import anndata as ad

from sc_browser.core.dataset import Dataset
from sc_browser.core.filter_state import FilterState


def _make_dataset():
    obs = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3", "c4"],
            "cluster": ["A", "A", "B", "B"],
            "condition": ["x", "y", "x", "y"],
            "sample": ["s1", "s1", "s2", "s2"],
        },
        index=["c1", "c2", "c3", "c4"],
    )

    var = pd.DataFrame(index=["g1", "g2"])
    X = np.arange(obs.shape[0] * var.shape[0]).reshape(obs.shape[0], var.shape[0])

    adata = ad.AnnData(X=X, obs=obs, var=var)

    ds = Dataset(
        name="TestDataset",
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


def test_subset_for_state_filters_correctly():

    ds = _make_dataset()

    state = FilterState(
        dataset_name=ds.name,
        view_id="subset",
        clusters=["A"],
        conditions=["x"],
        samples=[],
        genes=[],
    )

    sub = ds.subset_for_state(state)

    assert list(sub.adata.obs_names) == ["c1"]
    assert sub.adata.n_obs == 1
    assert list(sub.clusters) == ["A"]
    assert list(sub.conditions) == ["x"]


def test_sub_for_state_filters_samples():

    ds = _make_dataset()

    state = FilterState(
        dataset_name=ds.name,
        view_id="subset",
        clusters=[],
        conditions=[],
        samples=["s2"],
        genes=[],
    )

    sub = ds.subset_for_state(state)

    assert list(sub.adata.obs_names) == ["c3", "c4"]
    assert list(sub.clusters) == ["B", "B"]
