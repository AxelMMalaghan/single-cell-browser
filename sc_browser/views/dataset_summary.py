from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.graph_objs import Figure

from sc_browser.core.dataset import Dataset
from sc_browser.core.state import FilterState
from sc_browser.core.base_view import BaseView

class DatasetSummary(BaseView):
    """
    Dataset summary / inspector view

    - Respects current filters (clusters/conditions)
    - Shows bar charts of clister sizes and condition trees
    - Annotates n_cells, n_genes, and obs schema info
    """

    id = "summary"
    label = "Dataset summary"

    def compute_data(self, state: FilterState) -> pd.DataFrame:

        ds: Dataset = self.dataset.subset(
            clusters=self.clusters or None,
            conditions=self.conditions or None,
        )

        adata = ds.adata

        if adata.n_obs == 0:
            return pd.DataFrame(columns=["category", "count", "group_type"])

        cluster_counts = (
            ds.clusters.value_counts()
            .rename_axis("category")
            .reset_index(name="count")
        )
        cluster_counts["group_type"] = "Cluster"

        condition_counts = (
            ds.conditions.value_counts()
            .rename_axis("category")
            .reset_index(name="count")
        )

        condition_counts["group_type"] = "Condition"

        counts_df = pd.concat([cluster_counts, condition_counts], ignore_index=True)

        counts_df.attrs["n_cells"] = int(adata.n_obs)
        counts_df.attrs["n_genes"] = int(adata.n_vars)
        counts_df.attrs["dataset_name"] = ds.name

        obs = adata.obs
        schema_rows = []
        for col in obs.columns:
            series = obs[col]
            schema_rows.append(
                {
                    "column": str(col),
                    "dtype": str(series.dtype),
                    "n_unique": int(series.nunique()),
                }
            )
        schema_df = pd.DataFrame(schema_rows).sort_values("column")
        counts_df.attrs["obs_schema"] = schema_df

        return counts_df


    def render_figure(self, data: Any, state: FilterState) -> Figure:

        if data is None or data.empty:
            fig = go.Figure()
            fig.update_layout(
                title="No cells available for current filters",
                xaxis={"visible": False},
                yaxis={"visible": False},
            )
            return fig

        n_cells = data.attrs.get("n_cells", None)
        n_genes = data.attrs.get("n_genes", None)
        dataset_name = data.attrs.get("dataset_name", self.dataset.name)
        obs_schema: pd.DataFrame = data.attrs.get("obs_schema")

        # --- Main bar chart: cluster + condition sizes, faceted ---
        fig = px.bar(
            data,
            x="category",
            y="count",
            color="group_type",
            facet_col="group_type",
            facet_col_spacing=0.08,
            labels={"category": "Label", "count": "Cells", "group_type": ""},
        )

        fig.update_xaxes(tickangle=45)
        fig.update_layout(
            height=650,
            margin=dict(l=40, r=40, t=80, b=100),
            showlegend=False,
        )

        # --- Title / subtitle ---
        title_text = f"Dataset summary – {dataset_name}"
        subtitle_parts = []
        if n_cells is not None:
            subtitle_parts.append(f"Cells: {n_cells:,}")
        if n_genes is not None:
            subtitle_parts.append(f"Genes: {n_genes:,}")
        subtitle = " | ".join(subtitle_parts) if subtitle_parts else ""

        fig.update_layout(
            title={
                "text": f"{title_text}<br><sup>{subtitle}</sup>",
                "x": 0.0,
                "xanchor": "left",
            }
        )

        # --- Optional: embed a compact obs schema as text annotation ---
        if isinstance(obs_schema, pd.DataFrame) and not obs_schema.empty:
            # Limit to first 10 rows so it doesn't get insane
            head = obs_schema.head(10)
            lines = [
                f"{row['column']}  ({row['dtype']}, {row['n_unique']} levels)"
                for _, row in head.iterrows()
            ]
            if len(obs_schema) > len(head):
                lines.append(f"... (+{len(obs_schema) - len(head)} more columns)")

            schema_text = "obs columns:\n" + "\n".join(lines)

            fig.add_annotation(
                text=f"<b>Schema</b><br><span style='font-size:11px; white-space:pre'>{schema_text}</span>",
                xref="paper",
                yref="paper",
                x=1.02,
                y=1.0,
                xanchor="left",
                yanchor="top",
                showarrow=False,
                align="left",
            )

        return fig



