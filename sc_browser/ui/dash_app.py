# sc_browser/ui/dash_app.py (only the layout part changed)

from __future__ import annotations

from dash import Dash, dcc, html, Input, Output
import dash_bootstrap_components as dbc

from sc_browser.config import load_datasets
from sc_browser.core.state import FilterState
from sc_browser.core.view_registry import ViewRegistry
from sc_browser.views import ClusterView, ExpressionView, FeatureCountView, Dotplot, HeatmapView


def create_dash_app() -> Dash:
    # --- Load datasets via core config loader ----
    global_config, datasets = load_datasets("config/datasets.json")
    dataset_by_name = {ds.name: ds for ds in datasets}
    default_dataset = datasets[0]

    # --- Setup view registry ---
    registry = ViewRegistry()
    registry.register(ClusterView)
    registry.register(ExpressionView)
    registry.register(FeatureCountView)
    registry.register(Dotplot)
    registry.register(HeatmapView)

    app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

    # ------------- UI PIECES --------------

    navbar = dbc.Navbar(
        dbc.Container(
            fluid=True,
            children=[

                html.Div(
                    [
                        html.H2("Single-Cell Browser", className="mb-0"),
                        html.Small("Demo dataset explorer", className="text-muted"),
                    ],
                    className="d-flex flex-column justify-content-center",
                ),


                html.Div(
                    dcc.Dropdown(
                        id="dataset-select",
                        options=[{"label": ds.name, "value": ds.name} for ds in datasets],
                        value=default_dataset.name,
                        clearable=False,
                        placeholder="Select dataset",
                        style={"minWidth": "260px"},
                    ),
                    className="ms-auto d-flex align-items-center",
                    style={"width": "300px"},
                ),
            ],
        ),
        color="white",
        dark=False,
        className="shadow-sm scb-navbar",
    )

    # Left filter panel
    filter_panel = dbc.Card(
        [
            dbc.CardHeader("Filters", className="fw-semibold"),
            dbc.CardBody(
                [
                    html.Label("Filter clusters", className="form-label"),
                    dcc.Dropdown(
                        id="cluster-select",
                        options=[
                            {"label": c, "value": c}
                            for c in sorted(default_dataset.clusters.unique())
                        ],
                        multi=True,
                        placeholder="All clusters",
                        className="mb-3",
                    ),
                    html.Label("Filter conditions", className="form-label"),
                    dcc.Dropdown(
                        id="condition-select",
                        options=[
                            {"label": c, "value": c}
                            for c in sorted(default_dataset.conditions.unique())
                        ],
                        multi=True,
                        placeholder="All conditions",
                        className="mb-3",
                    ),
                    html.Label("Gene(s)", className="form-label"),
                    dcc.Dropdown(
                        id="gene-select",
                        options=[
                            {"label": g, "value": g}
                            for g in sorted(default_dataset.adata.var_names)
                        ],
                        multi=True,
                        placeholder="Select gene(s)",
                        className="mb-3",
                    ),
                    html.Hr(),
                    dbc.Checklist(
                        id="options-checklist",
                        options=[
                            {
                                "label": " Split by condition",
                                "value": "split_by_condition",
                            },
                        ],
                        value=[],
                        switch=True,
                    ),
                ]
            ),
        ],
        className="scb-sidebar",
    )

    # Right panel: tabs + main graph
    plot_panel = dbc.Card(
        [
            dbc.CardHeader(
                dcc.Tabs(
                    id="view-tabs",
                    value=ClusterView.id,
                    children=[
                        dcc.Tab(label=cls.label, value=cls.id)
                        for cls in registry.all_classes()
                    ],
                    className="scb-tabs",
                ),
                className="p-0",
            ),
            dbc.CardBody(
                dcc.Loading(
                    id="main-graph-loading",
                    type="default",
                    children=dcc.Graph(
                        id="main-graph",
                        style={"height": "650px"},
                        config={"responsive": True},
                    ),
                ),
                className="scb-main-body",
            ),
        ],
        className="scb-maincard",
    )

    app.layout = dbc.Container(
        fluid=True,
        className="scb-root",
        children=[
            navbar,
            dbc.Row(
                [
                    dbc.Col(filter_panel, md=3, className="mt-3"),
                    dbc.Col(plot_panel, md=9, className="mt-3"),
                ],
                className="gx-3",
            ),
        ],
    )

    # ---------- callbacks (your existing ones) ----------
    @app.callback(
        Output("cluster-select", "options"),
        Output("condition-select", "options"),
        Output("gene-select", "options"),
        Input("dataset-select", "value"),
    )
    def update_filters(dataset_name):
        ds = dataset_by_name[dataset_name]
        cluster_options = [{"label": c, "value": c} for c in sorted(ds.clusters.unique())]
        condition_options = [{"label": c, "value": c} for c in sorted(ds.conditions.unique())]
        gene_options = [{"label": g, "value": g} for g in sorted(ds.genes)]
        return cluster_options, condition_options, gene_options

    @app.callback(
        Output("main-graph", "figure"),
        Input("view-tabs", "value"),
        Input("dataset-select", "value"),
        Input("cluster-select", "value"),
        Input("condition-select", "value"),
        Input("gene-select", "value"),
        Input("options-checklist", "value"),
    )
    def update_main_graph(view_id, dataset_name, clusters, conditions, genes, options):
        ds = dataset_by_name[dataset_name]

        state = FilterState(
            genes=genes or [],
            clusters=clusters or [],
            conditions=conditions or [],
            split_by_condition="split_by_condition" in (options or []),
        )

        view = registry.create(view_id, ds)
        data = view.compute_data(state)
        fig = view.render_figure(data, state)
        return fig

    return app