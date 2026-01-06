from __future__ import annotations

import logging
from typing import List, Tuple

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from sc_browser.core.base_view import BaseView
from sc_browser.core.dataset import Dataset
from sc_browser.core.filter_profile import FilterProfile
from sc_browser.core.filter_state import FilterState
from sc_browser.services.differential_expression import DEConfig, run_de

logger = logging.getLogger(__name__)


class VolcanoPlotView(BaseView):
    """
    Volcano plot for differential expression.

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
        colour_scale=True,
    )

    def _choose_groups(
        self, state: FilterState, ds: Dataset
    ) -> Tuple[str, str, str | None]:
        groupby = ds.condition_key
        if groupby is None:
            raise ValueError(
                "No 'condition_key' configured for this dataset. Differential expression cannot run."
            )

        conditions = ds.conditions
        if conditions is None:
            raise ValueError("No condition values found in the dataset.")

        unique_conds: List[str] = sorted(conditions.astype(str).unique().tolist())
        if not unique_conds:
            raise ValueError("No condition values available for comparison.")

        if state.conditions:
            group1 = state.conditions[0]
            group2 = state.conditions[1] if len(state.conditions) > 1 else None
        else:
            group1 = unique_conds[0]
            group2 = None

        return groupby, group1, group2

    def compute_data(self, state: FilterState) -> pd.DataFrame:
        base_ds = self.dataset
        ds_de = base_ds.subset(
            clusters=state.clusters or None,
            samples=getattr(state, "samples", None) or None,
            cell_types=getattr(state, "cell_types", None) or None,
            conditions=None,
        )

        if ds_de.adata.n_obs == 0:
            return pd.DataFrame()

        try:
            groupby, group1, group2 = self._choose_groups(state, ds_de)
            config = DEConfig(
                dataset=ds_de, groupby=groupby, group1=group1, group2=group2
            )
            de_result = run_de(config)
            df = de_result.table.copy()
        except ValueError as e:
            # Handle expected configuration/data errors and log them
            logger.warning("VolcanoPlotView DE failed: %s", e)
            df = pd.DataFrame()
            df.attrs["error"] = str(e)
            return df
        except Exception:
            # Catch-all for unexpected engine failures
            logger.exception("Unexpected error in DE engine for VolcanoPlotView")
            return pd.DataFrame()

        if df.empty:
            return df

        # Cap -log10(p-value)
        p_col = "adj_pvalue" if "adj_pvalue" in df.columns else "pvalue"
        p_values = pd.Series(df[p_col]).astype(float)
        safe_p = p_values.replace(0, np.nan).to_numpy()
        neg_log_p = -np.log10(safe_p)
        neg_log_p_series = pd.Series(neg_log_p, index=p_values.index)
        max_finite = neg_log_p_series.replace([np.inf, -np.inf], np.nan).max()
        df["neg_log10_pvalue"] = neg_log_p_series.fillna(
            (max_finite or 50.0) + 1
        )

        # Thresholds and Significance
        df["significance"] = "Not Significant"
        sig_mask = p_values <= 0.05
        df.loc[sig_mask & (df["log2FC"] >= 1.0), "significance"] = "Upregulated"
        df.loc[sig_mask & (df["log2FC"] <= -1.0), "significance"] = "Downregulated"

        df.attrs["log_fc_threshold"] = 1.0
        df.attrs["pval_threshold"] = 0.05
        df.attrs["comparison"] = f"{group1} vs {'rest' if group2 is None else group2}"

        return df

    def render_figure(self, data: pd.DataFrame, state: FilterState) -> go.Figure:
        if data is None or data.empty:
            msg = (
                data.attrs.get("error", "No data to show")
                if data is not None
                else "No data to show"
            )
            return self.empty_figure(msg)

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

        fig.add_hline(
            y=-np.log10(pval_threshold),
            line_dash="dash",
            line_color="black",
            opacity=0.6,
        )
        fig.add_vline(
            x=log_fc_threshold, line_dash="dash", line_color="black", opacity=0.6
        )
        fig.add_vline(
            x=-log_fc_threshold, line_dash="dash", line_color="black", opacity=0.6
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
