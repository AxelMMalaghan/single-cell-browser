import anndata as ad
from pathlib import Path

path = Path("../data/FRC_LEC_IMM_FDC_TfH_GCB_sce_firstlook_500cells.h5ad")

adata = ad.read_h5ad(path)
adata.obs_names_make_unique()

print("obs columns:")
print(adata.obs.columns)

print("\nobsm keys:")
print(adata.obsm.keys())