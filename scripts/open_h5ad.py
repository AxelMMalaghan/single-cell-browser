import anndata

file_path = "../data/FRC_LEC_IMM_FDC_TfH_GCB_sce_firstlook_500cells.h5ad"

# Load the AnnData file (not in backed mode to inspect .obs)
adata = anndata.read_h5ad(file_path)

# Show available .obs columns
#print("obs columns:", list(adata.obs.columns))

# Preview first few rows of obs
#print("\nSample rows from .obs:")
#print(adata.obs.head())
#print(list(adata.obs.columns))
print(adata.obs.keys())