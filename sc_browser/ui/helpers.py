from __future__ import annotations

from typing import List, Tuple
from dash import html

from sc_browser.core.dataset import Dataset


def get_filter_dropdown_options(
    dataset: Dataset,
) -> Tuple[List[dict], List[dict], List[dict], List[dict], List[dict]]:
    cluster_options = [{"label": c, "value": c} for c in sorted(dataset.clusters.unique())]

    condition_series = dataset.conditions
    if condition_series is not None:
        condition_options = [{"label": c, "value": c} for c in sorted(condition_series.unique())]
    else:
        condition_options = []

    sample_options: List[dict] = []
    sample_key = dataset.obs_columns.get("sample")
    if sample_key and sample_key in dataset.adata.obs.columns:
        vals = dataset.adata.obs[sample_key].astype(str).unique()
        sample_options = [{"label": v, "value": v} for v in sorted(vals)]

    celltype_options: List[dict] = []
    celltype_key = dataset.obs_columns.get("cell_type")
    if celltype_key and celltype_key in dataset.adata.obs.columns:
        vals = dataset.adata.obs[celltype_key].astype(str).unique()
        celltype_options = [{"label": v, "value": v} for v in sorted(vals)]

    emb_options = [{"label": k, "value": k} for k in dataset.adata.obsm.keys()]

    return cluster_options, condition_options, sample_options, celltype_options, emb_options


def dataset_status(ds: Dataset) -> str:
    obs = ds.adata.obs
    obsm = ds.adata.obsm

    has_emb = ds.embedding_key in obsm
    has_cluster = ds.cluster_key in obs.columns
    has_condition = ds.condition_key in obs.columns

    if not has_emb:
        return f"❌ Missing embedding (obsm['{ds.embedding_key}'] not found)"
    if not has_cluster or not has_condition:
        return "⚠️ Missing or invalid cluster/condition mapping"
    return "✅ Ready"


def obs_preview_table(ds: Dataset, max_rows: int = 5) -> html.Table:
    obs = ds.adata.obs.copy()
    if obs.empty:
        return html.Table([html.Tr([html.Td("No obs columns")])])

    df = obs.reset_index().head(max_rows)
    columns = df.columns.tolist()
    rows = df.to_dict("records")

    header = html.Tr([html.Th(col) for col in columns])
    body = [html.Tr([html.Td(row[col]) for col in columns]) for row in rows]

    return html.Table(
        [html.Thead(header), html.Tbody(body)],
        className="table table-sm table-striped",
    )