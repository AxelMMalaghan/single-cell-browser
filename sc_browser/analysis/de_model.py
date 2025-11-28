from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import pandas as pd

from sc_browser.core.dataset import Dataset

@dataclass(frozen=True)
class DEConfig:
    """
    Immutable configuration for differential expression run.

    Keeps DE concerns separate from UI / Dash and from raw AnnData Cells
    """
    dataset: Dataset
    groupby: str
    group1: str
    group2: Optional[str] = None
    method: str = "wilcoxon"


@dataclass
class DEResult:
    """
    Canonical DE result: one row per gene

    Columns:
    - gene              gene-identifier
    - log2FC            log2 fold change
    - pvalue            raw p-value
    - adj_pvalue        multiple testing corrected p-value
    - comparison        human-readable label for comparison
    """

    config: DEConfig
    table: pd.DataFrame

    @property
    def signifcant(self) -> pd.DataFrame:
        """ Convenience: adjusted p < 0.05"""
        return self.table[self.table["adj_pvalue"] < 0.05].copy()

    def head(self, n: int = 10) -> pd.DataFrame:
        """Shortcut for quick inspection"""
        return self.table.head(n)