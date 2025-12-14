# Developer Documentation

This document is for developers working on the **Single‑Cell Browser** in an in‑house / research environment. It focuses on architecture, local setup, testing, and deployment hygiene rather than end‑user usage.

---

## 1. Project Overview

The Single‑Cell Browser is a Dash‑based web application for exploring single‑cell datasets (AnnData / Scanpy ecosystem). It provides:

* Dataset loading and validation
* Cached, filter‑driven subsetting
* Multiple views (cluster, expression, heatmap, volcano, etc.)
* Session‑scoped figure metadata and export
* Plotly‑based interactive and static outputs

The system is designed to be **state‑driven**, **pure where possible**, and **robust against malformed datasets**.

---

## 2. High‑Level Architecture

```
sc_browser/
├── core/           # Core domain logic (Dataset, FilterState)
├── views/          # Visualization logic (one class per view)
├── export/         # Figure rendering + export services
├── metadata_io/    # Session / figure metadata models
├── config/         # Dataset + UI configuration loading
├── ui/             # Dash layout + callbacks
└── app.py          # App entrypoint
```

### Design rules

* `core` must not depend on UI
* `views` depend on `core`, never the other way around
* `export` depends on `views` + `core`, not UI
* `config` is **input only** – avoid importing it into core logic

Violating these usually results in circular imports or brittle tests.

---

## 3. Core Data Model

### Dataset

`Dataset` is the central abstraction wrapping an `AnnData` object.

Responsibilities:

* Normalise and cache `.obs` columns (cluster, condition, sample, cell type)
* Efficient subsetting with caching
* Expression matrix extraction (sparse‑aware)
* Embedding access with validation

Key invariants:

* Datasets are **never mutated** after construction
* Subsets are cheap views backed by cached masks
* Missing embeddings raise explicit errors

### FilterState

`FilterState` is a **pure, serialisable** representation of user intent:

* Dataset + view context
* Selection filters (genes, clusters, conditions, etc.)
* Display options (dimensions, colour scale, 3D toggle)

Rules:

* Views must treat `FilterState` as immutable
* Anything persisted must round‑trip via `to_dict()` / `from_dict()`

---

## 4. Views

Each view class:

* Is responsible for **exactly one visualization type**
* Accepts a `Dataset` at construction
* Implements `render_figure(state, metadata)`

Views should:

* Call `dataset.subset_for_state(state)` as the *only* filtering entrypoint
* Fail gracefully on empty subsets
* Never access Dash components directly

Heavy computation belongs in the view layer, not callbacks.

---

## 5. Export Pipeline

Export flow:

```
FigureMetadata → ExportService → View.render_figure → Plotly Figure → File
```

Key points:

* `FigureMetadata.filter_state` may be a dict or `FilterState`
* ExportService must normalise this before rendering
* Plotly image export uses Kaleido

Exports must be deterministic and independent of Dash runtime state.

---

## 6. Metadata & Sessions

### FigureMetadata

Represents a *single saved figure*:

* Dataset key
* View ID
* FilterState
* View parameters
* Optional label + filename stem

### SessionMetadata

Represents a reporting session:

* Versioned schema
* Dataset config hash
* Ordered list of figures

⚠️ Important:

* `session_to_dict()` serialises nested dataclasses
* `session_from_dict()` **must reconstruct FilterState via `from_dict()`**

---

## 7. Local Development Setup

### Python

Recommended:

* Python 3.11+
* Virtual environment (`venv`, Poetry, or Conda)

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Run locally

```bash
python app.py
```

Dash runs on `http://localhost:8050` by default.

---

## 8. Testing

Tests are **first‑class** and required for deployment.

### Test categories

* Core logic (Dataset, FilterState)
* View computation (no Dash)
* Export + metadata round‑trips

### Run tests

```bash
pytest
```

Rules:

* Tests must not depend on Dash callbacks
* Use minimal in‑memory AnnData where possible
* No filesystem writes unless explicitly testing export

---

## 9. Docker & Deployment

The project includes a production‑ready `Dockerfile`.

### Build

```bash
docker build -t sc-browser .
```

### Run

```bash
docker run -p 8050:8050 sc-browser
```

Docker is the **reference deployment path** for in‑house use.

---

## 10. Coding Standards

* Black‑formatted Python
* Ruff for linting
* Explicit imports (avoid `__init__.py` side‑effects)
* No implicit global state

If you add a new view:

1. Add the class
2. Add at least one unit test
3. Register it in the view registry

---

## 11. Common Failure Modes

* Circular imports → usually caused by `__init__.py` imports
* Treating `FilterState` like a dict
* Forgetting to normalise obs column names
* Rendering large datasets without downsampling

If something feels "mysteriously broken", check imports first.
