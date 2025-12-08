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
    filter_profile = FilterProfile(
        clusters=True,
        conditions=True,
        samples=True,
        cell_types=True,
        genes=False,
        embedding=False,
        split_by_condition=False,
        is_3d=False,
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
        # For DE we must NOT apply the condition filter before running DE,
        # otherwise we can easily drop one of the comparison groups.
        base_ds = self.dataset

        # Apply only non-condition filters (hits subset cache)
        ds_de = base_ds.subset(
            clusters=state.clusters or None,
            samples=getattr(state, "samples", None) or None,
            cell_types=(
                getattr(state, "cell_types", None)
                or getattr(state, "celltypes", None)
                or None
            ),
            conditions=None,  # <-- crucial: keep all condition groups for DE
        )

        if ds_de.adata.n_obs == 0:
            return pd.DataFrame()

        # Pick groupby / group1 / group2 based on current state *and* DE dataset
        try:
            groupby, group1, group2 = self._choose_groups(state, ds_de)
        except ValueError:
            # No valid grouping → no DE to show
            return pd.DataFrame()

        config = DEConfig(
            dataset=ds_de,
            groupby=groupby,
            group1=group1,
            group2=group2,
            method="wilcoxon",
        )

        # Run DE; if groups are still invalid for some reason, fail gracefully
        try:
            de_result = run_de(config)
        except ValueError:
            return pd.DataFrame()

        df = de_result.table.copy()
        if df.empty:
            return df

        # ------------------------------------------------------------------
        # Derive volcano metrics and significance labels
        # ------------------------------------------------------------------
        # thresholds (could be later made configurable via UI)
        log_fc_threshold = 1.0
        pval_threshold = 0.05

        # Use adjusted p-value if available, else fall back to raw p-value
        if "adj_pvalue" in df.columns:
            p = df["adj_pvalue"].astype(float).copy()
        else:
            p = df["pvalue"].astype(float).copy()

        # Avoid log10(0) -> -inf
        p = p.replace(0, np.nan)
        df["neg_log10_pvalue"] = -np.log10(p)

        # Default: not significant
        df["significance"] = "Not Significant"

        # Up / down regulated masks
        sig_mask = p <= pval_threshold
        up_mask = (df["log2FC"] >= log_fc_threshold) & sig_mask
        down_mask = (df["log2FC"] <= -log_fc_threshold) & sig_mask

        df.loc[up_mask, "significance"] = "Upregulated"
        df.loc[down_mask, "significance"] = "Downregulated"

        # Attach metadata for render_figure
        df.attrs["log_fc_threshold"] = log_fc_threshold
        df.attrs["pval_threshold"] = pval_threshold
        if group2 is None:
            comparison = f"{group1} vs rest"
        else:
            comparison = f"{group1} vs {group2}"
        df.attrs["comparison"] = comparison

        return df

    def render_figure(self, data: pd.DataFrame, state: FilterState) -> go.Figure:

        if data is None or data.empty:
            return self.empty_figure("No data to show")

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
