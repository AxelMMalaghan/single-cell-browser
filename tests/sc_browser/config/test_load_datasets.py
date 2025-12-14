import json
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import pytest

from sc_browser.core.dataset import Dataset
from sc_browser.core.filter_state import FilterState
from sc_browser.config.loader import load_datasets


def _make_tiny_h5ad(tmp_path: Path) -> Path:
    obs = pd.DataFrame(
        {
            "cell_id": ["c1", "c2"],
            "cluster": ["A", "B"],
            "condition": ["x", "y"],
        },
        index=["c1", "c2"],
    )
    var = pd.DataFrame(index=["g1", "g2"])
    X = np.array([[1, 2], [3, 4]])

    adata = ad.AnnData(X=X, obs=obs, var=var)
    h5_path = tmp_path / "tiny.h5ad"
    adata.write_h5ad(h5_path)
    return h5_path


def test_load_datasets_from_config_dir(tmp_path):
    # Arrange: build config dir:
    # root/
    #   global.json
    #   datasets/
    #     dataset_0.json
    config_root = tmp_path / "config"
    datasets_dir = config_root / "datasets"
    datasets_dir.mkdir(parents=True)

    h5_path = _make_tiny_h5ad(tmp_path)

    # global.json
    global_json = {
        "ui_title": "Test Browser",
        "default_group": "Example",
        # data_root not needed for this test
    }
    (config_root / "global.json").write_text(json.dumps(global_json))

    # dataset config
    dataset_entry = {
        "name": "TinyDataset",
        "group": "ExampleGroup",
        "path": str(h5_path),
        "obs_columns": {
            "cell_id": "cell_id",
            "cluster": "cluster",
            "condition": "condition",
        },
        "embedding_key": "X_umap",
    }
    (datasets_dir / "dataset_0.json").write_text(json.dumps(dataset_entry))

    # Act
    global_config, datasets = load_datasets(config_root)

    # Assert: GlobalConfig wiring
    assert global_config.ui_title == "Test Browser"
    assert global_config.default_group == "Example"
    assert len(global_config.datasets) == 1

    # Dataset instances
    assert len(datasets) == 1
    ds = datasets[0]

    assert ds.name == "TinyDataset"
    assert ds.group == "ExampleGroup"
    assert ds.adata.n_obs == 2
    assert ds.adata.n_vars == 2

    # obs columns mapped correctly
    assert list(ds.clusters) == ["A", "B"]
    assert list(ds.conditions) == ["x", "y"]



def test_non_unique_obs_names_can_be_made_unique():
    X = np.random.rand(3, 2)
    adata = ad.AnnData(X=X)
    adata.obs_names = ["cell", "cell", "cell"]

    adata.obs_names_make_unique()

    assert len(set(adata.obs_names)) == adata.n_obs


def test_missing_embedding_key_detected():
    X = np.random.rand(5, 3)
    adata = ad.AnnData(X=X)
    adata.obs["cluster"] = ["A"] * 5
    adata.obs["condition"] = ["C"] * 5
    # No adata.obsm["X_umap"]

    ds = Dataset(
        name="test",
        group="g",
        adata=adata,
        cluster_key="cluster",
        condition_key="condition",
        embedding_key="X_umap",
        obs_columns=None,
    )

    with pytest.raises(ValueError, match="Embedding 'X_umap' not found"):
        ds.get_embedding()


def test_empty_subset_does_not_crash():
    X = np.random.rand(5, 3)
    adata = ad.AnnData(X=X)
    adata.obs["cluster"] = ["A"] * 5
    adata.obs["condition"] = ["C"] * 5
    adata.obsm["X_umap"] = np.random.rand(5, 2)

    ds = Dataset(
        name="test",
        group="g",
        adata=adata,
        cluster_key="cluster",
        condition_key="condition",
        embedding_key="X_umap",
        obs_columns=None,
    )

    state = FilterState(dataset_name="test", view_id="cluster", clusters=["DOES_NOT_EXIST"])
    sub = ds.subset_for_state(state)

    assert sub.adata.n_obs == 0

    # Should not crash to request embedding even if empty (embedding exists in obsm)
    emb = sub.get_embedding()
    assert emb.shape[0] == 0