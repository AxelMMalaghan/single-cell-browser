from __future__ import annotations

from typing import List

import pandas as pd
import scipy


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