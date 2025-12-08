from __future__ import annotations

from typing import Optional, Sequence

import scanpy as sc

from sc_browser.analysis.de_model import DEConfig, DEResult


def _validate_config(config: DEConfig) -> None:
    adata = config.dataset.adata
    groupby = config.groupby

    # Check grouping column exists
    if groupby not in adata.obs.columns:
        raise ValueError(
            f"groupby='{groupby}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    # Get values from that column and normalise to string
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


def _maybe_subset_genes(adata, config: DEConfig):
    """
    Optional: restrict DE to a subset of genes if config exposes a 'genes' field.

    This is purely an optimisation / convenience; if the config has no such field,
    we return the original AnnData unchanged.
    """
    genes: Optional[Sequence[str]] = getattr(config, "genes", None)
    if not genes:
        return adata

    gene_mask = adata.var_names.isin(genes)
    if not gene_mask.any():
        # No overlap: let scanpy handle the empty case later (will likely error),
        # or you could raise here if you want stricter behaviour.
        return adata

    # Return a view; caller decides if they copy
    return adata[:, gene_mask]


def run_de(config: DEConfig) -> DEResult:
    """
    Run a differential expression analysis using scanpy.rank_genes_groups.

    Notes:
    - If config.group1 is not None and config.group2 is None: group1 v. rest
    - If config.group2 is set: group1 v. group2 only (subset)
    """
    _validate_config(config)

    base_adata = config.dataset.adata
    groupby = config.groupby
    group1 = config.group1
    group2: Optional[str] = config.group2
    method = config.method

    # Optional gene restriction (if config has .genes)
    adata = _maybe_subset_genes(base_adata, config)

    # Group 1 v. rest
    if group2 is None:
        # Work on a copy to avoid mutating the shared Dataset AnnData
        work = adata.copy()
        sc.tl.rank_genes_groups(
            work,
            groupby=groupby,
            groups=[group1],
            reference="rest",
            method=method,
            use_raw=False,
        )
        comparison_label = f"{group1} v. rest"

    else:
        # Group 1 v. Group 2
        # First subset to only the two groups of interest
        mask = adata.obs[groupby].astype(str).isin([group1, group2])
        work = adata[mask].copy()

        # Normalise categories to str to match group1/group2
        work.obs[groupby] = work.obs[groupby].astype(str)

        g1 = str(group1)
        g2 = str(group2)

        sc.tl.rank_genes_groups(
            work,
            groupby=groupby,
            groups=[g1],
            reference=g2,
            method=method,
            use_raw=False,
        )
        comparison_label = f"{g1} v. {g2}"

    # Convert scanpy result to a flat DataFrame
    df = sc.get.rank_genes_groups_df(
        work,
        group=group1,
    )

    # Normalise schema
    df = df.rename(
        columns={
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
        raise RuntimeError(f"Missing columns from DE result: {missing}")

    df["comparison"] = comparison_label

    return DEResult(config=config, table=df[required + ["comparison"]])
