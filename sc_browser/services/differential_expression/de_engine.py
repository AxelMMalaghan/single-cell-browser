from __future__ import annotations

import logging
from functools import lru_cache
from typing import Literal, Optional, Sequence, cast

import numpy as np
import pandas as pd

from sc_browser.services.differential_expression.de_model import DEConfig, DEResult

logger = logging.getLogger(__name__)


def _validate_config(config: DEConfig) -> None:
    """Validate that groupby and group selections exist in the data."""
    adata = config.dataset.adata
    groupby = config.groupby

    if groupby not in adata.obs.columns:
        raise ValueError(
            f"groupby='{groupby}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    groups = adata.obs[groupby].astype(str).unique()

    if str(config.group1) not in groups:
        raise ValueError(
            f"group1='{config.group1}' not present in obs['{groupby}']. "
            f"Available groups: {list(groups)}"
        )

    if config.group2 is not None and str(config.group2) not in groups:
        raise ValueError(
            f"group2='{config.group2}' not present in obs['{groupby}']. "
            f"Available groups: {list(groups)}"
        )


def _maybe_subset_genes(adata, config: DEConfig):
    """Restrict DE to a subset of genes if configured."""
    genes: Optional[Sequence[str]] = getattr(config, "genes", None)
    if not genes:
        return adata

    gene_mask = adata.var_names.isin(genes)
    if not gene_mask.any():
        return adata

    return adata[:, gene_mask]


@lru_cache(maxsize=32)
def run_de(config: DEConfig) -> DEResult:
    """
    Run differential expression analysis using scanpy.rank_genes_groups.

    Includes heuristic check to prevent redundant log-transformation.
    """
    import scanpy as sc

    _validate_config(config)

    base_adata = config.dataset.adata
    groupby = config.groupby
    group1 = str(config.group1)
    group2: Optional[str] = str(config.group2) if config.group2 is not None else None
    method = cast(
        Literal["logreg", "t-test", "wilcoxon", "t-test_overestim_var"],
        config.method,
    )

    adata = _maybe_subset_genes(base_adata, config)

    # Prepare comparison subset
    if group2 is None:
        work = adata.copy()
        comparison_label = f"{group1} v. rest"
        ref = "rest"
    else:
        mask = adata.obs[groupby].astype(str).isin([group1, group2])
        work = adata[mask].copy()
        work.obs[groupby] = work.obs[groupby].astype(str)
        comparison_label = f"{group1} v. {group2}"
        ref = group2

    # FIX: Use a heuristic to detect if data is already log-transformed
    # even if 'log1p' is missing from .uns
    already_logged = "log1p" in work.uns
    if not already_logged:
        # Check max value in the first 100 cells/genes as a quick heuristic
        # If max value is > 50, it's likely raw counts.
        try:
            X = work.X
            if X is None:
                sample_max = 0
            else:
                sample_max = float(np.asarray(X[:100, :100]).max())
            if sample_max < 50:
                logger.info(
                    "Data appears already log-transformed (max < 50). Skipping sc.pp.log1p."
                )
                already_logged = True
        except Exception:
            pass

    if not already_logged:
        logger.info(f"Logarithmizing data for comparison: {comparison_label}")
        sc.pp.log1p(work)

    # Execute DE engine
    sc.tl.rank_genes_groups(
        work,
        groupby=groupby,
        groups=[group1],
        reference=ref,
        method=method,
        use_raw=False,
    )

    # Extract results
    df = sc.get.rank_genes_groups_df(work, group=group1)

    # Map to internal schema
    df = df.rename(
        columns={
            "names": "gene",
            "logfoldchanges": "log2FC",
            "pvals": "pvalue",
            "pvals_adj": "adj_pvalue",
        }
    )

    required = ["gene", "log2FC", "pvalue", "adj_pvalue"]
    df["comparison"] = comparison_label

    table = cast(pd.DataFrame, df[required + ["comparison"]])
    return DEResult(config=config, table=table)
