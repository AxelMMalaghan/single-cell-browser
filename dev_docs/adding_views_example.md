# Example: Adding a New View

This section walks through adding a **new visualization view** end‑to‑end (code + registration + tests). The goal is to keep the view:

* **Pure and deterministic** (same inputs → same outputs)
* **Independent of Dash** (no callbacks, no components)
* **Robust to empty subsets**

### 1.1 Choose a view concept

A view is one “unit of visualization” (e.g. cluster scatter, expression violin, volcano plot). For this example we’ll add:

* **Cell count by cluster** (simple but representative)

It will:

* subset the dataset using `FilterState`
* compute counts per cluster
* render a Plotly bar chart

### 1.2 Create the view class

Create a new file:

* `sc_browser/views/cluster_count_view.py`

```python
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import pandas as pd
import plotly.graph_objs as go

from sc_browser.core.dataset import Dataset
from sc_browser.core.filter_state import FilterState
from sc_browser.export.model import FigureMetadata


class ClusterCountView:
    """Bar chart showing number of cells per cluster for the current filter state."""

    # Convention (match your existing view patterns)
    id: str = "cluster_count"
    label: str = "Cluster count"

    def __init__(self, dataset: Dataset) -> None:
        self._dataset = dataset

    def render_figure(self, state: FilterState, metadata: FigureMetadata) -> go.Figure:
        """Render a Plotly figure for the given state."""
        sub = self._dataset.subset_for_state(state)

        # If no cluster series exists, render an explicit message
        clusters = sub.clusters
        if clusters is None:
            fig = go.Figure()
            fig.add_annotation(
                text="No cluster mapping configured for this dataset.",
                showarrow=False,
            )
            fig.update_xaxes(visible=False)
            fig.update_yaxes(visible=False)
            return fig

        # Empty subset → empty chart with message
        if sub.adata.n_obs == 0:
            fig = go.Figure()
            fig.add_annotation(
                text="No cells match the current filters.",
                showarrow=False,
            )
            fig.update_xaxes(visible=False)
            fig.update_yaxes(visible=False)
            return fig

        # Compute counts per cluster
        counts = (
            clusters.value_counts(dropna=False)
            .rename_axis("cluster")
            .reset_index(name="n_cells")
            .sort_values("n_cells", ascending=False)
        )

        fig = go.Figure(
            data=[
                go.Bar(
                    x=counts["cluster"].astype(str),
                    y=counts["n_cells"].astype(int),
                )
            ]
        )

        # Minimal layout: consistent title and axes
        title = metadata.label or "Cells per cluster"
        fig.update_layout(title=title, xaxis_title="Cluster", yaxis_title="Cells")

        return fig
```

Notes:

* The view uses **only** `Dataset` + `FilterState` + `FigureMetadata`.
* It handles missing cluster mapping and empty subsets explicitly.

### 1.3 Register the view

Where you register views (typically a registry module), add the new class.

Example (adjust to your actual registry API):

```python
# sc_browser/core/view_registry.py (or equivalent)

from sc_browser.views.cluster_count_view import ClusterCountView

registry.register(ClusterCountView)
```

If your registry expects a mapping, do:

```python
VIEWS = {
    "cluster_count": ClusterCountView,
    # ... other views
}
```

### 1.4 Add a unit test

Create:

* `tests/sc_browser/views/test_cluster_count_view.py`

```python
from __future__ import annotations

import numpy as np
import anndata as ad

from sc_browser.core.dataset import Dataset
from sc_browser.core.filter_state import FilterState
from sc_browser.export.model import FigureMetadata
from sc_browser.views.cluster_count_view import ClusterCountView


def test_cluster_count_view_renders_counts():
    X = np.random.rand(5, 3)
    adata = ad.AnnData(X=X)
    adata.obs["cluster"] = ["A", "A", "B", "B", "B"]

    ds = Dataset(
        name="test",
        group="g",
        adata=adata,
        cluster_key="cluster",
        condition_key=None,
        embedding_key=None,
        obs_columns=None,
    )

    view = ClusterCountView(ds)

    state = FilterState(dataset_name="test", view_id="cluster_count")

    md = FigureMetadata.from_runtime(
        figure_id="fig-0001",
        dataset_key="test",
        view_id="cluster_count",
        state=state,
        view_params={},
        label="Counts",
        file_stem="fig-0001",
    )

    fig = view.render_figure(state, md)

    # Basic assertion: 2 bars (A and B)
    assert len(fig.data) == 1
    assert list(fig.data[0].x) == ["B", "A"]  # sorted by count desc
    assert list(fig.data[0].y) == [3, 2]
```

This test:

* builds a tiny in‑memory AnnData
* constructs a Dataset and view
* asserts the computed Plotly trace matches expected values

### 1.5 Wire into the UI

The UI typically uses view IDs + labels to populate dropdowns.

Once the registry sees the new view class, it should appear automatically if your UI:

* calls `registry.all_classes()` and builds options from `.label` and `.id`

If you maintain an explicit allow‑list per dataset/profile, add `"cluster_count"` there too.

### 1.6 Checklist

Before merging:

* [ ] View renders on empty subsets
* [ ] View renders when required mappings are missing
* [ ] Unit test added
* [ ] Registered in view registry
* [ ] No UI callback changes needed unless you added new controls

---