from __future__ import annotations

from dash import Dash, dcc, html, Input, Output
import dash_bootstrap_components as dbc

from sc_browser.config import load_datasets
from sc_browser.core.state import FilterState
from sc_browser.core.view_registry import ViewRegistry
from sc_browser.views import ClusterView

def create_dash_app() -> Dash:

    # --- Load datasets via core config loader ----
    global_config, datasets = load_datasets("config/datasets.json")
    dataset_by_name = {ds.name: ds for ds in datasets}
    default_dataset = datasets[0]

    # --- Setup view registry ---
    registry = ViewRegistry()
    registry.register(ClusterView)

    app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

    app.layout = dbc.Container(
        fluid=True,
        children=[
            dbc.Row([
                dbc.Col(
                    width=3,
                    children=[
                        html.H3("Single-cell Browser"),
                        dcc.Dropdown(
                            id="dataset-select",
                            options=[{"label": ds.name, "value": ds.name} for ds in datasets],
                            value=default_dataset.name,
                            clearable=False,
                            style={"margin-bottom": "1rem"},
                        ),
                        dcc.Dropdown(
                            id="cluster-select",
                            options=[
                                {"label": c, "value": c}
                                for c in sorted(default_dataset.clusters.unique())
                            ],
                            multi=True,
                            placeholder="Filter clusters (optional)",
                            style={"margin-bottom": "1rem"},
                        ),
                        dcc.Dropdown(
                            id="condition-select",
                            options=[
                                {"label": c, "value": c}
                                for c in sorted(default_dataset.conditions.unique())
                            ],
                            multi=True,
                            placeholder="Filter conditions (optional)",
                            style={"margin-bottom": "1rem"},
                        ),
                        dbc.Checklist(
                            id="options-checklist",
                            options=[
                                {"label": "Split by condition", "value": "split_by_condition"},
                            ],
                            value=[],
                            style={"margin-bottom": "1rem"},
                        ),
                    ],
                ),
                dbc.Col(
                    width=9,
                    children=[
                        dcc.Tabs(
                            id="view-tabs",
                            value=ClusterView.id,
                            children=[
                                dcc.Tab(label=cls.label, value=cls.id)
                                for cls in registry.all_classes()
                            ],
                        ),
                        dcc.Graph(
                            id="main-graph",
                            style={"height": "600px"},  # or 500 / 700 / whatever you like
                            config={"responsive": True},
                        ),
                    ],
                ),
            ])
        ],
    )

    # --- Callbacks ---

    @app.callback(
        Output("cluster-select", "options"),
        Output("condition-select", "options"),
        Input("dataset-select", "value")
    )

    def update_filters(dataset_name):

        ds = dataset_by_name[dataset_name]
        cluster_options = [{"label": c, "value": c} for c in sorted(ds.clusters.unique())]
        condition_options = [{"label": c, "value": c} for c in sorted(ds.conditions.unique())]
        return cluster_options, condition_options

    @app.callback(
        Output("main-graph", "figure"),
        Input("view-tabs", "value"),
        Input("dataset-select", "value"),
        Input("cluster-select", "value"),
        Input("condition-select", "value"),
        Input("options-checklist", "value"),
    )

    def update_main_graph(view_id, dataset_name, clusters, conditions, options):
        ds = dataset_by_name[dataset_name]

        state = FilterState(clusters=clusters or [],
                            conditions=conditions or [],
                            split_by_condition="split_by_condition" in (options or []),
        )

        view = registry.create(view_id, ds)
        data = view.compute_data(state)
        fig = view.render_figure(data, state)
        return fig

    return app
