import numpy as np
import pandas as pd
import anndata as ad

from sc_browser.core.dataset import Dataset


def _make_dataset_for_expr():
    obs = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3"],
            "cluster": ["A", "A", "B"],
            "condition": ["x", "y", "x"],
            "sample": ["s1", "s1", "s2"],
            "cell_type": ["T", "T", "B"],
        },
        index=["c1", "c2", "c3"],
    )
    var = pd.DataFrame(index=["g1", "g2", "g3"])
    X = np.array(
        [
            [1, 0, 5],  # c1
            [0, 2, 6],  # c2
            [3, 4, 0],  # c3
        ]
    )

    adata = ad.AnnData(X=X, obs=obs, var=var)

    ds = Dataset(
        name="ExprDataset",
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
            "cell_type": "cell_type",
        },
        file_path=None,
    )
    return ds


def test_expression_matrix_filters_genes_and_caches():
    ds = _make_dataset_for_expr()

    # First call
    genes = ["g1", "g3"]
    df1 = ds.expression_matrix(genes)

    assert list(df1.columns) == ["g1", "g3"]
    assert list(df1.index) == ["c1", "c2", "c3"]
    assert df1.loc["c1", "g1"] == 1
    assert df1.loc["c2", "g3"] == 6

    # Second call with same genes should hit cache (same object)
    df2 = ds.expression_matrix(genes)
    assert df1 is df2

    # Request gene that doesn't exist -> empty df (but indexed by obs_names)
    df3 = ds.expression_matrix(["missing_gene"])
    assert df3.empty
    assert list(df3.index) == ["c1", "c2", "c3"]


def test_subset_cache_returns_same_object_for_same_filters():
    ds = _make_dataset_for_expr()

    sub1 = ds.subset(
        clusters=["A"],
        conditions=["x"],
        samples=["s1"],
        cell_types=["T"],
    )
    sub2 = ds.subset(
        clusters=["A"],
        conditions=["x"],
        samples=["s1"],
        cell_types=["T"],
    )

    # Cached subset → identity
    assert sub1 is sub2

    # Different filters → different Dataset instance
    sub3 = ds.subset(clusters=["B"])
    assert sub3 is not sub1
    assert list(sub3.adata.obs_names) == ["c3"]