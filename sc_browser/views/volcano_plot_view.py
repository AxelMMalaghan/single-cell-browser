from __future__ import annotations

from typing import Any, List, Tuple

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from sc_browser.analysis import DEConfig, run_de
from sc_browser.core.base_view import BaseView
from sc_browser.core.filter_state import FilterState, FilterProfile
from sc_browser.core.dataset import Dataset


class VolcanoPlotView(BaseView):
    """
    Volcano plot for differential expression analysis.

    X: log2 fold change
    Y: -log10(p-value)
    Colour: up / down / not significant
    """

    id = "volcano"
    label = "Volcano"
    filter_profile = FilterProfile(clusters=True,
                                   conditions=False,
                                   samples=True,
                                   cell_types=True,
                                   genes=False
    )

    def _choose_groups(self, state: FilterState, ds: Dataset) -> Tuple[str, str, str | None]:
        """
        Decide which obs column and groups to compare.

        MVP behaviour:
        - use condition_key as groupby
        - if user has picked ≥1 conditions: first = group1
        - if user has picked ≥2 conditions: second = group2
        - else: group1 = first condition in dataset, group2 = None (vs rest)
        """
        groupby = ds.condition_key
        if groupby is None:
            raise ValueError("No condition_key configured for DE")

        conditions = ds.conditions
        if conditions is None:
            raise ValueError("No condition values available for DE")

        unique_conds: List[str] = sorted(conditions.astype(str).unique().tolist())
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
        # Apply filters (hits subset cache); DE runs on the filtered dataset
        ds: Dataset = self.dataset.subset(
            clusters=state.clusters or None,
            conditions=state.conditions or None,
            samples=getattr(state, "samples", None) or None,
            cell_types=getattr(state, "cell_types", None) or None,
        )

        if ds.adata.n_obs == 0:
            return pd.DataFrame()

        # Pick groupby / group1 / group2 based on current state *and* filtered dataset
        try:
            groupby, group1, group2 = self._choose_groups(state, ds)
        except ValueError:
            # No valid grouping → no DE to show
            return pd.DataFrame()

        config = DEConfig(
            dataset=ds,          # use subset, not full dataset
            groupby=groupby,
            group1=group1,
            group2=group2,
            method="wilcoxon",
        )

        de_result = run_de(config)
        df = de_result.table.copy()

        if df.empty:
            return df

        # Pre-compute -log10(p) and significance flags
        eps = 1e-300
        p = df["pvalue"].to_numpy()
        p_safe = np.clip(p, eps, 1.0)
        df["neg_log10_pvalue"] = -np.log10(p_safe)

        # Thresholds – could be made configurable later
        log_fc_threshold = 1.0
        pval_threshold = 0.05

        df["significance"] = "Not Significant"
        up = (df["log2FC"] > log_fc_threshold) & (df["pvalue"] < pval_threshold)
        down = (df["log2FC"] < -log_fc_threshold) & (df["pvalue"] < pval_threshold)

        df.loc[up, "significance"] = "Upregulated"
        df.loc[down, "significance"] = "Downregulated"

        # Store thresholds and comparison meta for renderer
        df.attrs["log_fc_threshold"] = log_fc_threshold
        df.attrs["pval_threshold"] = pval_threshold

        # be defensive: table may not have comparison column if run_de changes
        comparison_col = "comparison"
        if comparison_col in df.columns and not df[comparison_col].empty:
            df.attrs["comparison"] = df[comparison_col].iloc[0]
        else:
            df.attrs["comparison"] = f"{group1} vs {group2 or 'rest'}"

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
