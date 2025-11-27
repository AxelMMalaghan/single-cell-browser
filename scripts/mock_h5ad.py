import numpy as np
import pandas as pd
import anndata as ad
from pathlib import Path

n_cells = 100
n_genes = 50

X = np.random.poisson(lam=1.0, size=(n_cells, n_genes))

obs = pd.DataFrame(
    {
        "cluster":   np.random.choice(["C0", "C1", "C2"], size=n_cells),
        "condition": np.random.choice(["ctrl", "stim"], size=n_cells),
    },
    index=[f"cell_{i}" for i in range(n_cells)],
)

var = pd.DataFrame(index=[f"gene_{j}" for j in range(n_genes)])

umap = np.random.normal(size=(n_cells, 2))

adata = ad.AnnData(X=X, obs=obs, var=var)
adata.obsm["X_umap"] = umap

Path("data").mkdir(exist_ok=True)
adata.write_h5ad("data/demo.h5ad")
print("wrote data/demo.h5ad", adata.shape)