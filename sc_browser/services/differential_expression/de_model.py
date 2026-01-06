from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Tuple, cast

import pandas as pd

from sc_browser.core.dataset import Dataset


@dataclass(frozen=True)
class DEConfig:
    """Immutable configuration for a differential expression run."""

    dataset: Dataset
    groupby: str
    group1: str
    group2: Optional[str] = None
    method: str = "wilcoxon"
    genes: Optional[Tuple[str, ...]] = None

    def with_genes(self, genes: Sequence[str]) -> "DEConfig":
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
        """Thresholds: Adjusted P-value <= 0.05 and |log2FC| >= 1.0."""
        if self.table.empty:
            return self.table
        mask = (self.table["adj_pvalue"] <= 0.05) & (self.table["log2FC"].abs() >= 1.0)
        return cast(pd.DataFrame, self.table[mask])

    def head(self, n: int = 10) -> pd.DataFrame:
        return self.table.head(n)
