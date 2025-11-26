import numpy as np
import pandas as pd
import anndata as ad
from pathlib import Path

# Tiny toy dataset
n_cells = 10000
n_genes = 5000

X = np.random.poisson(1.0, size=(n_cells, n_genes)).astype(np.float32)

obs = pd.DataFrame(
    {
        "X__cluster": np.random.choice(["C0", "C1", "C2"], size=n_cells),
        "condition": np.random.choice(["Ctrl", "Stim"], size=n_cells),
    },
    index=[f"cell{i}" for i in range(n_cells)],
)

var = pd.DataFrame(
    index=[f"gene{i}" for i in range(n_genes)]
)

# Fake UMAP embedding
obsm = {
    "X_umap": np.random.normal(size=(n_cells, 2)).astype(np.float32)
}

adata = ad.AnnData(X=X, obs=obs, var=var, obsm=obsm)

out_path = Path("../data/demo.h5ad")
out_path.parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(out_path)

print("Wrote", out_path.resolve())
