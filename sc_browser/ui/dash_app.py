from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple, Optional
import base64
import logging

import dash
from dash import Dash, dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
import pandas as pd
import plotly.graph_objs as go

from sc_browser.config.io import load_datasets
from sc_browser.config.write import save_dataset_config
from sc_browser.core.state import FilterState
from sc_browser.core.view_registry import ViewRegistry
from sc_browser.views import (
    ClusterView,
    ExpressionView,
    FeatureCountView,
    DotplotView,
    HeatmapView,
    VolcanoPlotView,
    DatasetSummary,
)

logger = logging.getLogger(__name__)


def get_filter_dropdown_options(dataset) -> Tuple[List[dict], List[dict]]:
    """
    Return options for cluster and condition dropdowns.

    NOTE: genes are handled via a separate server-side search callback to avoid
    sending 30k+ options to the browser.
    """
    cluster_options = [{"label": c, "value": c} for c in sorted(dataset.clusters.unique())]
    condition_options = [{"label": c, "value": c} for c in sorted(dataset.conditions.unique())]
    return cluster_options, condition_options


def _build_view_registry() -> ViewRegistry:
    registry = ViewRegistry()
    registry.register(ClusterView)
    registry.register(ExpressionView)
    registry.register(FeatureCountView)
    registry.register(DotplotView)
    registry.register(HeatmapView)
    registry.register(VolcanoPlotView)
    registry.register(DatasetSummary)
    return registry


def _build_navbar(datasets, global_config) -> dbc.Navbar:
    title = getattr(global_config, "ui_title", "scB++")
    subtitle = getattr(global_config, "subtitle", "Interactive Dataset Explorer")
    return dbc.Navbar(
        dbc.Container(
            fluid=True,
            children=[
                html.Div(
                    [
                        html.H2(title, className="mb-0"),
                        html.Small(subtitle, className="text-muted", id="navbar-subtitle"),
                    ],
                    className="d-flex flex-column justify-content-center",
                ),
                html.Div(
                    dcc.Dropdown(
                        id="dataset-select",
                        options=[{"label": ds.name, "value": ds.name} for ds in datasets],
                        value=datasets[0].name if datasets else None,
                        clearable=False,
                        placeholder="Select dataset",
                        style={"minWidth": "260px"},
                        className="scb-dataset-dropdown",
                    ),
                    className="ms-auto d-flex align-items-center",
                    style={"width": "300px"},
                ),
            ],
        ),
        color="light",
        dark=False,
        className="shadow-sm scb-navbar",
    )


def _build_filter_panel(default_dataset) -> dbc.Card:
    cluster_options, condition_options = get_filter_dropdown_options(default_dataset)
    return dbc.Card(
        [
            dbc.CardHeader("Filters", className="fw-semibold"),
            dbc.CardBody(
                [
                    html.Label("Filter clusters", className="form-label"),
                    dcc.Dropdown(
                        id="cluster-select",
                        options=cluster_options,
                        multi=True,
                        placeholder="All clusters",
                        className="mb-3",
                    ),
                    html.Label("Filter conditions", className="form-label"),
                    dcc.Dropdown(
                        id="condition-select",
                        options=condition_options,
                        multi=True,
                        placeholder="All conditions",
                        className="mb-3",
                    ),
                    html.Label("Gene(s)", className="form-label"),
                    dcc.Dropdown(
                        id="gene-select",
                        options=[],  # populated via server-side search callback
                        multi=True,
                        placeholder="Type to search genes",
                        className="mb-3",
                    ),
                    html.Hr(),
                    dbc.Checklist(
                        id="options-checklist",
                        options=[{"label": " Split by condition", "value": "split_by_condition"}],
                        value=[],
                        switch=True,
                    ),
                    html.Hr(),
                    html.Div(
                        [
                            html.Div(
                                default_dataset.name,
                                id="sidebar-dataset-name",
                                className="fw-semibold",
                            ),
                            html.Div(
                                [
                                    html.Span(
                                        f"{default_dataset.adata.n_obs} cells · {default_dataset.adata.n_vars} genes",
                                        id="sidebar-dataset-meta",
                                        className="text-muted",
                                    ),
                                ],
                                className="small",
                            ),
                        ],
                        className="scb-dataset-summary mt-1",
                    ),
                ]
            ),
        ],
        className="scb-sidebar",
    )


def _build_plot_panel(registry: ViewRegistry) -> dbc.Card:
    return dbc.Card(
        [
            dbc.CardHeader(
                dcc.Tabs(
                    id="view-tabs",
                    value=ClusterView.id,
                    children=[
                        dcc.Tab(label=cls.label, value=cls.id) for cls in registry.all_classes()
                    ],
                    className="scb-tabs",
                ),
                className="p-0",
            ),
            dbc.CardBody(
                [
                    dcc.Loading(
                        id="main-graph-loading",
                        type="default",
                        children=dcc.Graph(
                            id="main-graph",
                            style={"height": "650px"},
                            config={"responsive": True},
                        ),
                    ),
                    html.Div(
                        [
                            dbc.Button(
                                "Download data (CSV)",
                                id="download-data-btn",
                                color="secondary",
                                size="sm",
                                className="mt-2 me-2",
                            ),
                            dcc.Download(id="download-data"),
                        ],
                        className="d-flex justify-content-end",
                    ),
                ],
                className="scb-main-body",
            ),
        ],
        className="scb-maincard",
    )


def _dataset_status(ds) -> str:
    """Return a human-readable status string for the dataset."""
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


def _obs_preview_table(ds, max_rows: int = 5) -> html.Table:
    """Build a small HTML table preview of .obs."""
    obs = ds.adata.obs.copy()
    if obs.empty:
        return html.Table([html.Tr([html.Td("No obs columns")])])

    df = obs.reset_index().head(max_rows)
    columns = df.columns.tolist()
    rows = df.to_dict("records")

    header = html.Tr([html.Th(col) for col in columns])
    body = [
        html.Tr([html.Td(row[col]) for col in columns]) for row in rows
    ]

    return html.Table(
        [html.Thead(header), html.Tbody(body)],
        className="table table-sm table-striped",
    )


def _build_dataset_manager_panel() -> dbc.Container:
    """
    Datasets screen: manage / inspect / map datasets + import.

    Uses the global dataset-select as the source of truth.
    Layout:
      - Top row: status/import (left), mapping form (right)
      - Bottom row: full-width obs preview
    """
    status_card = dbc.Card(
        [
            dbc.CardHeader("Current dataset"),
            dbc.CardBody(
                [
                    html.Div(
                        id="dm-current-dataset",
                        className="mb-1 fw-semibold",
                    ),
                    html.Div(
                        id="dm-status-text",
                        className="mb-1",
                    ),
                    html.Div(
                        id="dm-summary-text",
                        className="text-muted mb-3",
                    ),
                    html.Hr(),
                    html.H6("Import dataset (.h5ad)", className="mt-1"),
                    dcc.Upload(
                        id="dm-upload",
                        children=html.Div(
                            [
                                "Drag and drop or ",
                                html.A("select a .h5ad file"),
                            ]
                        ),
                        multiple=False,
                        className="scb-upload border rounded p-2 text-center",
                    ),
                    html.Div(
                        id="dm-import-status",
                        className="small text-muted mt-2",
                    ),
                ]
            ),
        ],
        className="h-100",
    )

    mapping_card = dbc.Card(
        [
            dbc.CardHeader("Mapping"),
            dbc.CardBody(
                [
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    html.Label("Cluster column", className="form-label"),
                                    dcc.Dropdown(id="dm-cluster-key", className="mb-2"),
                                ],
                                md=6,
                            ),
                            dbc.Col(
                                [
                                    html.Label("Condition column", className="form-label"),
                                    dcc.Dropdown(id="dm-condition-key", className="mb-2"),
                                ],
                                md=6,
                            ),
                        ],
                        className="gx-3",
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    html.Label("Sample column", className="form-label"),
                                    dcc.Dropdown(id="dm-sample-key", className="mb-2"),
                                ],
                                md=6,
                            ),
                            dbc.Col(
                                [
                                    html.Label("Cell type column", className="form-label"),
                                    dcc.Dropdown(id="dm-celltype-key", className="mb-2"),
                                ],
                                md=6,
                            ),
                        ],
                        className="gx-3",
                    ),
                    html.Label("Embedding key", className="form-label mt-2"),
                    dcc.Dropdown(id="dm-embedding-key", className="mb-3"),
                    dbc.Button(
                        "Save mapping",
                        id="dm-save-btn",
                        color="primary",
                        size="sm",
                        className="mt-1",
                    ),
                    html.Span(
                        id="dm-save-status",
                        className="ms-2 small text-muted",
                    ),
                ]
            ),
        ],
        className="h-100",
    )

    obs_card = dbc.Card(
        [
            dbc.CardHeader("Obs preview"),
            dbc.CardBody(
                html.Div(id="dm-obs-preview"),
                className="p-2",
            ),
        ],
        className="mt-3",
    )

    return dbc.Container(
        fluid=True,
        children=[
            dbc.Row(
                [
                    dbc.Col(status_card, md=4, className="mt-3"),
                    dbc.Col(mapping_card, md=8, className="mt-3"),
                ],
                className="gx-3",
            ),
            dbc.Row(
                [
                    dbc.Col(obs_card, md=12),
                ],
                className="gx-3",
            ),
        ],
        className="scb-datasets-view",
    )


def create_dash_app() -> Dash:
    config_root = Path("config")

    global_config, datasets = load_datasets(config_root)
    if not datasets:
        raise RuntimeError("No datasets were loaded from config")

    # These will be mutated by the import callback via `nonlocal`
    dataset_by_name: Dict[str, object] = {ds.name: ds for ds in datasets}

    registry = _build_view_registry()

    app = Dash(__name__, external_stylesheets=[dbc.themes.FLATLY])

    navbar = _build_navbar(datasets, global_config)
    filter_panel = _build_filter_panel(datasets[0])
    plot_panel = _build_plot_panel(registry)
    dataset_manager = _build_dataset_manager_panel()

    app.layout = dbc.Container(
        fluid=True,
        className="scb-root",
        children=[
            navbar,
            dcc.Tabs(
                id="page-tabs",
                value="explore",
                children=[
                    dcc.Tab(
                        label="Explore",
                        value="explore",
                        children=[
                            dbc.Row(
                                [
                                    dbc.Col(filter_panel, md=3, className="mt-3"),
                                    dbc.Col(plot_panel, md=9, className="mt-3"),
                                ],
                                className="gx-3",
                            ),
                        ],
                    ),
                    dcc.Tab(
                        label="Datasets",
                        value="datasets",
                        children=[
                            dataset_manager,
                        ],
                    ),
                ],
                className="mt-2",
            ),
        ],
    )

    # -------------------------------------------------------------------------
    # Explore callbacks
    # -------------------------------------------------------------------------
    @app.callback(
        Output("sidebar-dataset-name", "children"),
        Output("sidebar-dataset-meta", "children"),
        Input("dataset-select", "value"),
    )
    def update_sidebar_dataset_summary(dataset_name: str):
        ds = dataset_by_name[dataset_name]
        name = ds.name
        meta = f"{ds.adata.n_obs} cells · {ds.adata.n_vars} genes"
        return name, meta



    @app.callback(
        Output("cluster-select", "value"),
        Output("condition-select", "value"),
        Output("gene-select", "value"),
        Input("dataset-select", "value"),
    )
    def reset_filters_on_dataset_change(dataset_name: str):
        # When the dataset changes, clear all filters.
        return [], [], []


    @app.callback(
        Output("cluster-select", "options"),
        Output("condition-select", "options"),
        Input("dataset-select", "value"),
    )
    def update_filters(dataset_name: str):
        return get_filter_dropdown_options(dataset_by_name[dataset_name])

    @app.callback(
        Output("gene-select", "options"),
        Input("gene-select", "search_value"),
        State("dataset-select", "value"),
        State("gene-select", "value"),
    )
    def update_gene_options(search_value, dataset_name, selected_genes):
        ds = dataset_by_name[dataset_name]
        all_genes = list(ds.genes)
        gene_set = set(all_genes)

        selected_genes = selected_genes or []

        # Remove stale genes from previous dataset
        selected_genes = [g for g in selected_genes if g in gene_set]

        # No search performed yet → only return selected genes
        if not search_value or not search_value.strip():
            return [{"label": g, "value": g} for g in selected_genes]

        query = search_value.lower()
        matches = [g for g in all_genes if query in g.lower()]
        suggestions = matches[:100]

        union = list(dict.fromkeys(selected_genes + suggestions))
        return [{"label": g, "value": g} for g in union]

    @app.callback(
        Output("main-graph", "figure"),
        Input("view-tabs", "value"),
        Input("dataset-select", "value"),
        Input("cluster-select", "value"),
        Input("condition-select", "value"),
        Input("gene-select", "value"),
        Input("options-checklist", "value"),
    )
    def update_main_graph(

        view_id: str,
        dataset_name: str,
        clusters: List[str] | None,
        conditions: List[str] | None,
        genes: List[str] | None,
        options: List[str] | None,
    ):
        ds = dataset_by_name[dataset_name]


        if genes:
            gene_set = set(ds.genes)
            genes = [g for g in genes if g in gene_set]

        state = FilterState(
            genes=genes or [],
            clusters=clusters or [],
            conditions=conditions or [],
            split_by_condition="split_by_condition" in (options or []),
        )

        view = registry.create(view_id, ds)


        try:
            data = view.compute_data(state)
            fig = view.render_figure(data, state)
            return fig
        except ValueError as e:
            logger.warning("View error (%s): %s", view_id, e)
            msg = str(e)
        except Exception:
            logger.exception("Unexpected error in view %s", view_id)
            msg = "An unexpected error occurred while rendering this view."

        label = getattr(view, "label", view_id)
        fallback = go.Figure()
        fallback.add_annotation(
            text=msg,
            showarrow=False,
            xref="paper",
            yref="paper",
            x=0.5,
            y=0.5,
            font=dict(size=14),
        )
        fallback.update_layout(
            title=f"{label} – unavailable",
            xaxis={"visible": False},
            yaxis={"visible": False},
            margin=dict(l=40, r=40, t=60, b=40),
        )
        return fallback

    @app.callback(
        Output("download-data", "data"),
        Input("download-data-btn", "n_clicks"),
        State("view-tabs", "value"),
        State("dataset-select", "value"),
        State("cluster-select", "value"),
        State("condition-select", "value"),
        State("gene-select", "value"),
        State("options-checklist", "value"),
        prevent_initial_call=True,
    )
    def download_current_data(
        n_clicks,
        view_id: str,
        dataset_name: str,
        clusters,
        conditions,
        genes,
        options,
    ):
        from dash import dcc as dash_dcc  # avoid name clash with global dcc

        ds = dataset_by_name[dataset_name]
        state = FilterState(
            genes=genes or [],
            clusters=clusters or [],
            conditions=conditions or [],
            split_by_condition="split_by_condition" in (options or []),
        )
        view = registry.create(view_id, ds)
        data = view.compute_data(state)
        if not isinstance(data, pd.DataFrame) or data.empty:
            return None
        filename = f"{view_id}_{dataset_name.replace(' ', '_')}.csv"
        return dash_dcc.send_data_frame(data.to_csv, filename, index=False)

    # -------------------------------------------------------------------------
    # Dataset manager callbacks (all driven by dataset-select)
    # -------------------------------------------------------------------------

    @app.callback(
        Output("dm-status-text", "children"),
        Output("dm-summary-text", "children"),
        Output("dm-cluster-key", "options"),
        Output("dm-cluster-key", "value"),
        Output("dm-condition-key", "options"),
        Output("dm-condition-key", "value"),
        Output("dm-sample-key", "options"),
        Output("dm-sample-key", "value"),
        Output("dm-celltype-key", "options"),
        Output("dm-celltype-key", "value"),
        Output("dm-embedding-key", "options"),
        Output("dm-embedding-key", "value"),
        Output("dm-obs-preview", "children"),
        Output("dm-current-dataset", "children"),
        Input("dataset-select", "value"),
    )
    def update_dataset_manager(dataset_name: str):
        ds = dataset_by_name[dataset_name]

        status = _dataset_status(ds)
        summary = f"{ds.adata.n_obs} cells · {ds.adata.n_vars} genes"

        obs_cols = list(ds.adata.obs.columns)
        obs_options = [{"label": c, "value": c} for c in obs_cols]

        cluster_val = ds.cluster_key if ds.cluster_key in ds.adata.obs.columns else None
        condition_val = ds.condition_key if ds.condition_key in ds.adata.obs.columns else None
        sample_val = ds.obs_columns.get("sample")
        celltype_val = ds.obs_columns.get("cell_type")

        emb_keys = list(ds.adata.obsm.keys())
        emb_options = [{"label": k, "value": k} for k in emb_keys]
        emb_val = ds.embedding_key if ds.embedding_key in ds.adata.obsm else (emb_keys[0] if emb_keys else None)

        obs_preview = _obs_preview_table(ds)
        current_label = f"Current dataset: {ds.name}"

        return (
            status,
            summary,
            obs_options,
            cluster_val,
            obs_options,
            condition_val,
            obs_options,
            sample_val,
            obs_options,
            celltype_val,
            emb_options,
            emb_val,
            obs_preview,
            current_label,
        )

    @app.callback(
        Output("dm-save-status", "children"),
        Input("dm-save-btn", "n_clicks"),
        State("dataset-select", "value"),
        State("dm-cluster-key", "value"),
        State("dm-condition-key", "value"),
        State("dm-sample-key", "value"),
        State("dm-celltype-key", "value"),
        State("dm-embedding-key", "value"),
        prevent_initial_call=True,
    )
    def save_dataset_mapping(
        n_clicks: int,
        dataset_name: str,
        cluster_key: Optional[str],
        condition_key: Optional[str],
        sample_key: Optional[str],
        celltype_key: Optional[str],
        embedding_key: Optional[str],
    ):
        if not n_clicks:
            raise dash.exceptions.PreventUpdate

        ds = dataset_by_name[dataset_name]

        if cluster_key:
            ds.cluster_key = cluster_key
        if condition_key:
            ds.condition_key = condition_key
        if embedding_key:
            ds.embedding_key = embedding_key

        obs_cols = dict(ds.obs_columns)
        if sample_key:
            obs_cols["sample"] = sample_key
        if celltype_key:
            obs_cols["cell_type"] = celltype_key
        obs_cols["cluster"] = ds.cluster_key
        obs_cols["condition"] = ds.condition_key
        ds.obs_columns = obs_cols

        path = save_dataset_config(ds, config_root)
        return f"Saved mapping to {path.relative_to(config_root.parent)}"

    @app.callback(
        Output("dm-import-status", "children"),
        Output("dataset-select", "options"),
        Output("dataset-select", "value"),
        Input("dm-upload", "contents"),
        State("dm-upload", "filename"),
        State("dataset-select", "value"),
        prevent_initial_call=True,
    )
    def import_dataset(contents, filename, current_dataset_name):
        nonlocal global_config, datasets, dataset_by_name

        if not contents or not filename:
            raise dash.exceptions.PreventUpdate

        if not filename.endswith(".h5ad"):
            opts = [{"label": ds.name, "value": ds.name} for ds in datasets]
            current = current_dataset_name or (datasets[0].name if datasets else None)
            return (
                f"Unsupported file type for '{filename}'. Please upload a .h5ad file.",
                opts,
                current,
            )

        # Decode uploaded file
        try:
            content_type, content_string = contents.split(",", 1)
            decoded = base64.b64decode(content_string)
        except Exception as e:
            logger.exception("Failed to decode uploaded file")
            opts = [{"label": ds.name, "value": ds.name} for ds in datasets]
            current = current_dataset_name or (datasets[0].name if datasets else None)
            return (
                f"Failed to decode uploaded file: {e}",
                opts,
                current,
            )

        # Save to data/ directory next to config/
        data_dir = config_root.parent / "data"
        data_dir.mkdir(parents=True, exist_ok=True)
        out_path = data_dir / filename

        try:
            with out_path.open("wb") as f:
                f.write(decoded)
        except Exception as e:
            logger.exception("Failed to write uploaded file to disk")
            opts = [{"label": ds.name, "value": ds.name} for ds in datasets]
            current = current_dataset_name or (datasets[0].name if datasets else None)
            return (
                f"Failed to save file to {out_path}: {e}",
                opts,
                current,
            )

        # Reload all datasets
        try:
            global_config, new_datasets = load_datasets(config_root)
        except Exception as e:
            logger.exception("Reloading datasets after import failed")
            opts = [{"label": ds.name, "value": ds.name} for ds in datasets]
            current = current_dataset_name or (datasets[0].name if datasets else None)
            return (
                f"File saved to {out_path}, but reloading datasets failed: {e}",
                opts,
                current,
            )

        datasets = new_datasets
        dataset_by_name = {ds.name: ds for ds in datasets}

        # Try to match uploaded file to a Dataset
        imported_ds_name = None
        for ds in datasets:
            fp = getattr(ds, "file_path", None)
            if fp is not None and fp.name == filename:
                imported_ds_name = ds.name
                break

        if imported_ds_name is None and datasets:
            imported_ds_name = datasets[-1].name

        opts = [{"label": ds.name, "value": ds.name} for ds in datasets]
        selected = imported_ds_name or (datasets[0].name if datasets else None)

        status_msg = f"Imported '{filename}'"
        if imported_ds_name:
            status_msg += f" as dataset '{imported_ds_name}'."
        else:
            status_msg += " (could not match to a specific dataset; using existing list.)"

        return (
            status_msg,
            opts,
            selected,
        )

    return app
