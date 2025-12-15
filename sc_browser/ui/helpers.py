from __future__ import annotations

from typing import List, Tuple
from dash import dash_table, html


from sc_browser.validation.dataset_validation import validate_dataset
from sc_browser.validation.errors import ValidationError
from sc_browser.core.dataset import Dataset

def get_filter_dropdown_options(
    dataset: Dataset,
) -> Tuple[List[dict], List[dict], List[dict], List[dict], List[dict]]:
    # --- Cluster options (defensive for datasets without cluster mapping) ---
    cluster_series = dataset.clusters
    if cluster_series is not None:
        cluster_options = [
            {"label": c, "value": c}
            for c in sorted(cluster_series.unique())
        ]
    else:
        cluster_options = []

    # --- Condition options (already defensive) ---
    condition_series = dataset.conditions
    if condition_series is not None:
        condition_options = [
            {"label": c, "value": c}
            for c in sorted(condition_series.unique())
        ]
    else:
        condition_options = []

    # --- Sample options (optional obs column) ---
    sample_options: List[dict] = []
    sample_key = dataset.obs_columns.get("sample")
    if sample_key and sample_key in dataset.adata.obs.columns:
        vals = dataset.adata.obs[sample_key].astype(str).unique()
        sample_options = [{"label": v, "value": v} for v in sorted(vals)]

    # --- Cell type options (optional obs column) ---
    celltype_options: List[dict] = []
    celltype_key = dataset.obs_columns.get("cell_type")
    if celltype_key and celltype_key in dataset.adata.obs.columns:
        vals = dataset.adata.obs[celltype_key].astype(str).unique()
        celltype_options = [{"label": v, "value": v} for v in sorted(vals)]

    # --- Embedding options (always from obsm keys) ---
    emb_options = [{"label": k, "value": k} for k in dataset.adata.obsm.keys()]

    return cluster_options, condition_options, sample_options, celltype_options, emb_options


def dataset_status(ds: Dataset):
    obs = ds.adata.obs
    obsm = ds.adata.obsm

    has_emb = ds.embedding_key in obsm
    has_cluster = ds.cluster_key in obs.columns if ds.cluster_key is not None else False
    has_condition = ds.condition_key in obs.columns if ds.condition_key is not None else False

    if not has_emb:
        return html.Span(
            [
                html.Strong("Status: "),
                "Missing embedding ",
                html.Code(f"obsm['{ds.embedding_key}']"),
                " not found.",
            ],
            className="dm-status dm-status-error",
        )

    if not has_cluster or not has_condition:
        missing_bits = []
        if not has_cluster:
            missing_bits.append("cluster")
        if not has_condition:
            missing_bits.append("condition")
        missing_text = "/".join(missing_bits)

        return html.Span(
            [
                html.Strong("Status: "),
                f"Missing or invalid {missing_text} mapping.",
            ],
            className="dm-status dm-status-warn",
        )

    return html.Span(
        [
            html.Strong("Status: "),
            "Ready.",
        ],
        className="dm-status dm-status-ok",
    )



def obs_preview_table(ds: Dataset, max_rows: int = 20):
    """
    Build a styled Dash DataTable showing a preview of .obs.
    """
    df = ds.adata.obs.copy()
    df = df.reset_index().head(max_rows)

    return dash_table.DataTable(
        data=df.to_dict("records"),
        columns=[{"name": c, "id": c} for c in df.columns],

        # ---- FONT + LOOK & FEEL ----
        style_table={
            "overflowX": "auto",
        },
        style_as_list_view=True,
        style_cell={
            "fontFamily": 'system-ui, -apple-system, BlinkMacSystemFont, "SF Pro Text", "Segoe UI", sans-serif',
            "fontSize": "12px",
            "padding": "6px 8px",
            "border": "none",
            "textAlign": "left",
            "minWidth": "80px",
            "maxWidth": "260px",
            "whiteSpace": "nowrap",
            "textOverflow": "ellipsis",
        },
        style_header={
            "fontFamily": 'system-ui, -apple-system, BlinkMacSystemFont, "SF Pro Text", "Segoe UI", sans-serif',
            "fontSize": "12px",
            "fontWeight": "600",
            "backgroundColor": "#f3f4f6",
            "borderBottom": "1px solid #e5e7eb",
        },
        style_data={
            "borderBottom": "1px solid #e5e7eb",
        },

        # optional: small tweaks
        page_size=max_rows,
        sort_action="native",
        filter_action="none",
    )




def warn_on_invalid_datasets(datasets: Iterable[Dataset], logger: logging.Logger) -> None:
    """
    Validate datasets and log warnings for any that are not deployment-ready.

    This is intentionally warn-only (in-house friendly): the app still runs,
    but you get actionable signals in logs immediately after load/import.
    """
    for ds in datasets:
        try:
            validate_dataset(ds)
        except ValidationError as e:
            logger.warning(
                "Dataset %r validation failed: %s",
                ds.name,
                "; ".join(f"{issue.code}: {issue.message}" for issue in e.issues),
            )
