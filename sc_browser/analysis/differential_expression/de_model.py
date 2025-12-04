from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple, Sequence

import pandas as pd

from sc_browser.core.dataset import Dataset


@dataclass(frozen=True)
class DEConfig:
    """
    Immutable configuration for a differential expression run.

    Keeps DE concerns separate from UI / Dash and from raw AnnData.
    Designed to be hashable / cache-friendly.
    """

    dataset: Dataset
    groupby: str
    group1: str
    group2: Optional[str] = None
    method: str = "wilcoxon"

    # Optional: restrict DE to this subset of genes
    # This is a performance lever – if None, DE runs on all genes.
    genes: Optional[Tuple[str, ...]] = None

    def with_genes(self, genes: Sequence[str]) -> "DEConfig":
        """
        Return a copy of this config with a restricted gene set.

        Useful if you want to build the base config once, then vary
        the gene subset (e.g. user-selected genes, HVGs, etc.).
        """
        # Normalise to a tuple of unique gene names in order of appearance
        unique_genes = tuple(dict.fromkeys(str(g) for g in genes))
        return DEConfig(
            dataset=self.dataset,
            groupby=self.groupby,
            group1=self.group1,
            group2=self.group2,
            method=self.method,
            genes=unique_genes,
        )


@dataclass
class DEResult:
    """
    Canonical DE result: one row per gene.

    Columns:
    - gene              gene identifier
    - log2FC            log2 fold change
    - pvalue            raw p-value
    - adj_pvalue        multiple testing corrected p-value
    - comparison        human-readable label for comparison
    """

    config: DEConfig
    table: pd.DataFrame

    @property
    def significant(self) -> pd.DataFrame:
        """Convenience: adjusted p < 0.05."""
        return self.table[self.table["adj_pvalue"] < 0.05].copy()

    # Backwards-compatible alias for the old typo, in case anything used it
    @property
    def signifcant(self) -> pd.DataFrame:  # type: ignore[override]
        return self.significant

    def head(self, n: int = 10) -> pd.DataFrame:
        """Shortcut for quick inspection."""
        return self.table.head(n)
