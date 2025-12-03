from __future__ import annotations

from typing import Any, List

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from sc_browser.analysis import DEConfig
from sc_browser.core.base_view import BaseView
from sc_browser.core.state import FilterState
from sc_browser.analysis import run_de


class VolcanoPlotView(BaseView):
    """
    Volcano plot for differential expression analysis

    X: log2 fold change
    Y: -log10(p-value)
    Colour: up / down / not significant
    """

    id = "volcano"
    label = "Volcano"

    def _choose_groups(self, state: FilterState) -> tuple[str, str | None, str]:
        """
        Decide which obs column and groups to compare

        For MVP:
        - use condition_key as groupby
        - if user has picked 1 or more conditions, first = group1
        - if user has picked 2 or more conditions, second = group2
        - else: group1 = first condition in dataset, group2 = None (vs rest)

        :param state: state of the FilterState - what the user has toggled for
        :return: groups
        """

        groupby = self.dataset.condition_key

        conditions = self.dataset.conditions.astype(str)
        unique_conds: List[str] = sorted(conditions.unique().tolist())

        if not unique_conds:
            raise ValueError("No condition values available for DE")

        if state.conditions:
            group1 = state.conditions[0]
            group2 = state.conditions[1] if len(state.conditions) > 1 else None
        else:
            group1 = unique_conds[0]
            group2 = None

        return groupby, group1, group2


    def compute_data(self, state: FilterState) -> pd.DataFrame:

        ds = self.dataset.subset(
            clusters=state.clusters or None,
            conditions=state.conditions or None,
        )

        # Pick groupby / group1 / group2 based on current state
        groupby, group1, group2 = self._choose_groups(state)

        config = DEConfig(
            dataset=ds,          # use subset, not full dataset
            groupby=groupby,
            group1=group1,
            group2=group2,
            method="wilcoxon",
        )

        de_result = run_de(config)
        df = de_result.table.copy()

        # Pre compute -log10(p) and significance flags
        # Avoid log10(0)
        eps = 1e-300
        p = df["pvalue"].to_numpy()
        p_safe = np.clip(p, eps, 1.0)
        df["neg_log10_pvalue"] = -np.log10(p_safe)

        log_fc_threshold = 1.0
        pval_threshold = 0.05

        df["significance"] = "Not Significant"
        up = (df["log2FC"] > log_fc_threshold) & (df["pvalue"] < pval_threshold)
        down = (df["log2FC"] < -log_fc_threshold) & (df["pvalue"] < pval_threshold)

        df.loc[up, "significance"] = "Upregulated"
        df.loc[down, "significance"] = "Downregulated"

        # Store thresholds for the renderer (so we don't hardcode twice)
        df.attrs["log_fc_threshold"] = log_fc_threshold
        df.attrs["pval_threshold"] = pval_threshold
        df.attrs["comparison"] = de_result.table["comparison"].iloc[0]

        return df

    def render_figure(self, data: pd.DataFrame, state: FilterState) -> Any:
        if data is None or data.empty:
            fig = go.Figure()
            fig.update_layout(
                title="No DE results available for current selection",
                xaxis={"visible": False},
                yaxis={"visible": False},
            )
            return fig

        log_fc_threshold = data.attrs.get("log_fc_threshold", 1.0)
        pval_threshold = data.attrs.get("pval_threshold", 0.05)
        comparison = data.attrs.get("comparison", "")

        fig = px.scatter(
            data,
            x="log2FC",
            y="neg_log10_pvalue",
            color="significance",
            hover_data={
                "gene": True,
                "log2FC": True,
                "pvalue": True,
                "adj_pvalue": True,
                "neg_log10_pvalue": True,
            },
            color_discrete_map={
                "Not Significant": "lightgray",
                "Upregulated": "red",
                "Downregulated": "blue",
            },
        )

        # Threshold lines
        fig.add_hline(
            y=-np.log10(pval_threshold),
            line_dash="dash",
            line_color="black",
            opacity=0.6,
        )
        fig.add_vline(
            x=log_fc_threshold,
            line_dash="dash",
            line_color="black",
            opacity=0.6,
        )
        fig.add_vline(
            x=-log_fc_threshold,
            line_dash="dash",
            line_color="black",
            opacity=0.6,
        )

        fig.update_layout(
            height=600,
            margin=dict(l=40, r=40, b=40, t=40),
            xaxis_title="log2(Fold Change)",
            yaxis_title="-log10(p-value)",
            legend_title="Significance",
            title=f"Volcano plot ({comparison})",
        )

        return fig

