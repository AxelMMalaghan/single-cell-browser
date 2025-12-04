# scripts/make_mock_dataset.py

from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad


def main() -> None:
    # project root = parent of this file's directory
    root = Path(__file__).resolve().parent.parent
    data_dir = root / "data"
    data_dir.mkdir(exist_ok=True)

    # ---- basic sizes ----
    n_cells = 1000000
    n_genes = 5000

    rng = np.random.default_rng(42)

    # ---- expression matrix X (cells x genes) ----
    # Poisson counts so feature/count plots look sane
    X = rng.poisson(size=(n_cells, n_genes)).astype("float32")

    # ---- obs: per-cell metadata ----
    obs = pd.DataFrame(
        {
            # this must match "cluster_key": "cluster" in demo1dataset.json
            "cluster": rng.choice(["C0", "C1", "C2", "C3"], size=n_cells),
            # this must match "condition_key": "condition" in demo1dataset.json
            "condition": rng.choice(["ctrl", "stim"], size=n_cells),
        },
        index=[f"cell_{i}" for i in range(n_cells)],
    )

    # ---- var: per-gene metadata ----
    var = pd.DataFrame(
        index=[f"gene_{j}" for j in range(n_genes)]
    )

    # ---- embedding: goes into obsm["X_umap"] ----
    # 2D UMAP-like coordinates
    umap = rng.normal(size=(n_cells, 2)).astype("float32")

    # ---- build AnnData ----
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obsm["X_umap"] = umap

    out_path = data_dir / "demo3.h5ad"
    adata.write_h5ad(out_path)

    print(f"Wrote: {out_path}")
    print("obs columns:", list(adata.obs.columns))
    print("obsm keys:", list(adata.obsm.keys()))
    print("shape:", adata.shape)


if __name__ == "__main__":
    main()