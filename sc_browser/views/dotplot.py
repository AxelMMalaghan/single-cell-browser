from __future__ import annotations

from typing import Any

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from sc_browser.core.dataset import Dataset
from sc_browser.core.view_base import BaseView
from sc_browser.core.state import FilterState

class Dotplot(BaseView):
    """
    Dot plot:
    - x-axis: genes
    - y-axis: cluster (or cluster-condition)
    - dot-size %:of cells expressing the gene in that group
    - dot-color: relative mean expression in that group
    """

    id = "dotplot"
    label = "Dot Plot"

    def compute(self, data: pd.DataFrame, state: FilterState) -> pd.DataFrame:
        adata = self.dataset.adata
        obs = self.obs.copy()

        # --- mask by cluster / condition