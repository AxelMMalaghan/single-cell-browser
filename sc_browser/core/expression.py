import numpy as np
import pandas as pd
import scipy

from typing import List
from anndata import AnnData


def get_dense_matrix(self, adata: AnnData, layer: str = None, max_cells: int = 10000) -> np.ndarray:
    """
    Safely returns a dense matrix from .X or a specified layer. Raises is matrix is too large to convert.
    :param self:
    :param adata:
    :param layer:
    :param max_cells:
    :return:
    """
    X = adata.layers[layer] if layer else adata.X

    # only warn or block if it's sparse
    if hasattr(X, "toarray"):
        nume1 = X.shape[0] * X.shape[1]
        if nume1 > max_cells * X.shape[1]:
            raise MemoryError("Attempted to densify {X.shape[0]}×{X.shape[1]} matrix "
                              f"({nume1:,} elements). This may crash your system.")
        return X.toarray()

    return X



def extract_expression_matrix(self, adata, genes: List[str]) -> pd.DataFrame:
    """
    Extracts expression matrix from gene expression matrix.
    :param self:
    :param adata:
    :param genes:
    :return:
    """

    X = adata[:, genes].X
    var_names = adata[:, genes].var_names

    if scipy.sparse.issparse(X):
        df = pd.DataFrame.sparse.from_spmatrix(X, columns=var_names)
    else:
        df = pd.DataFrame(X, columns=var_names)

    df.index = adata.obs_names
    return df


