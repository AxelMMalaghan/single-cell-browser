import os

import anndata as ad

# 1. Setup path
script_dir = os.path.dirname(os.path.abspath(__file__))
# Adjusting to go up one level from 'scripts/' to project root, then into 'data/'
path = os.path.join(script_dir, "../data/KH_combined_2023-Jan-11.h5ad")

print(f"Loading: {os.path.abspath(path)}")
adata = ad.read_h5ad(path)

# 2. FIX THE WARNING: Make variable (gene) names unique
if not adata.var_names.is_unique:
    print("\n⚠️  Fixing duplicate gene names...")
    adata.var_names_make_unique()

# 3. Inspect .obs (Cell Metadata)
print("\n=== Cell Metadata Columns (.obs) ===")
# This prints the LIST of columns available (e.g., 'cluster', 'sample', 'condition')
print(list(adata.obs.columns))

# Optional: Print the first few rows to see actual data
# print(adata.obs.head())

# 4. Inspect .obsm (Embeddings)
print("\n=== Embeddings Keys (.obsm) ===")
# This prints keys like 'X_umap', 'X_pca', etc.
print(list(adata.obsm.keys()))
