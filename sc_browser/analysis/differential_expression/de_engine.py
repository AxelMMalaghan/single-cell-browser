from __future__ import annotations

from typing import Optional

import scanpy as sc

from sc_browser.analysis.differential_expression.de_model import DEConfig, DEResult

def _validate_config(config: DEConfig) -> None:

    adata = config.dataset.adata
    groupby = config.groupby

    # Check grouping column exists
    if groupby not in adata.obs.columns:
        raise ValueError(
            f"groupby='{groupby}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    # Get **values** from that column and normalise to string
    groups = adata.obs[groupby].astype(str).unique()

    # Validate group1
    if config.group1 not in groups:
        raise ValueError(
            f"group1='{config.group1}' not present in obs['{groupby}']. "
            f"Available groups: {list(groups)}"
        )

    # Validate group2 only if it's actually set
    if config.group2 is not None and config.group2 not in groups:
        raise ValueError(
            f"group2='{config.group2}' not present in obs['{groupby}']. "
            f"Available groups: {list(groups)}"
        )


def run_de(config: DEConfig) -> DEResult:
    """
    Run a differential expression analysis using scanpy.rank_genes_groups.

    Notes:
    - If config.group1 is None: group1 v. rest
    - If config.group2 is set: group1 v. group2 only (subset)
    """
    _validate_config(config)

    adata = config.dataset.adata
    groupby = config.groupby
    group1 = config.group1
    group2: Optional[str] = config.group2
    method = config.method

    # Work on a copy to not mutate the AnnData
    if group2 is None:
        # Group 1 v. rest
        work = adata.copy()
        sc.tl.rank_genes_groups(
            work,
            groupby=groupby,
            groups=[group1],
            reference="rest",
            method=method,
            use_raw=False
        )
        comparison_label = f"{group1} v. rest"
    else:
        # Group 1 v. Group 2
        mask = adata.obs[groupby].isin([group1, group2])
        work = adata[mask].copy()
        work.obs[groupby] = work.obs[groupby].astype(str)

        group1 = str(config.group1)
        group2 = str(config.group2) if config.group2 is not None else None

        sc.tl.rank_genes_groups(
            work,
            groupby=groupby,
            groups=[group1],
            reference=group2,
            method=method,
            use_raw=False
        )
        comparison_label = f"{group1} v. {group2}"

    # Convert scanpy result to a flat DataFrame
    df = sc.get.rank_genes_groups_df(
        work,
        group=config.group1,
    )

    # Normalise schema //TODO: standardise schema?
    df = df.rename(columns={
        "names": "gene",
        "logfoldchanges": "log2FC",
        "pvals": "pvalue",
        "pvals_adj": "adj_pvalue",
        }
    )

    # Ensure required columns exist
    required = ["gene", "log2FC", "pvalue", "adj_pvalue"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise RuntimeError(f"Missing columns: {missing}")

    df["comparison"] = comparison_label

    return DEResult(config=config, table=df[required+ ["comparison"]])


