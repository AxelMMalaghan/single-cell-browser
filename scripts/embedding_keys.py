import anndata as ad



print("obsm keys:", adata.obsm.keys())
for key in adata.obsm.keys():
    print(key, adata.obsm[key].shape)