import anndata as ad
import numpy as np
import pandas as pd

from sc_browser.config.dataset_loader import from_config
from sc_browser.core.configs import DatasetConfig


def _make_h5ad_with_dupes(tmp_path):
    """Helper to create a small .h5ad with duplicate obs/var names."""
    obs = pd.DataFrame(
        {
            "cell_id": ["c1", "c1", "c2"],
            "cluster": ["A", "A", "B"],
            "condition": ["x", "y", "x"],
        },
        index=pd.Index(["c1", "c1", "c2"]),  # duplicate obs_names on purpose
    )

    var = pd.DataFrame(
        index=pd.Index(["g1", "g1", "g2"]),  # duplicate var_names on purpose
    )

    X = np.arange(obs.shape[0] * var.shape[0]).reshape(obs.shape[0], var.shape[0])

    adata = ad.AnnData(X=X, obs=obs, var=var)

    h5_path = tmp_path / "dupe_test.h5ad"
    adata.write_h5ad(h5_path)

    return h5_path


def test_from_config_normalises_obs_and_var_names(tmp_path):
    # Arrange: create test h5ad
    h5_path = _make_h5ad_with_dupes(tmp_path)

    raw = {
        "name": "TestDataset",
        "group": "Example",
        "path": str(h5_path),
        "obs_columns": {
            "cell_id": "cell_id",
            "cluster": "cluster",
            "condition": "condition",
        },
        "embedding_key": "X_umap",  # doesn't actually exist, but fine for wiring
    }

    cfg = DatasetConfig.from_raw(raw, source_path=h5_path, index=0)

    # Act
    ds = from_config(cfg)

    # Assert: obs/var names are unique and dimensions preserved
    adata = ds.adata

    assert adata.n_obs == 3
    assert adata.n_vars == 3

    assert adata.obs_names.is_unique
    assert adata.var_names.is_unique

    # Cluster/condition wiring
    assert ds.cluster_key == "cluster"
    assert ds.condition_key == "condition"
    assert ds.clusters is not None
    assert ds.conditions is not None
    assert list(ds.clusters.unique()) == ["A", "B"]
    assert set(ds.conditions.unique()) == {"x", "y"}

    # Dataset metadata
    assert ds.name == "TestDataset"
    assert ds.group == "Example"
    assert ds.file_path == h5_path
