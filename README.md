# Single-Cell Browser

An in-house, config-driven web application for interactive exploration of single-cell datasets (`.h5ad`), built on **Dash**, **Scanpy**, and **Plotly**.

The browser is designed to be:
- **Dataset-agnostic** (no hard-coded schema assumptions)
- **Robust to partial or imperfect metadata**
- **Configurable via JSON**, not code changes
- **Deployable on an internal server** (Docker + Gunicorn)

This is not a hosted SaaS product. It is intended for **controlled, internal use** in a research environment.

---

## High-Level Architecture

.h5ad files
↓
DatasetConfig (JSON)
↓
Dataset Loader
↓
Dataset abstraction
↓
FilterState
↓
Views (cluster, expression, DE, etc.)
↓
Dash UI (Plotly figures)

Key design principles:
- **All data access flows through the `Dataset` abstraction**
- **Views never touch AnnData directly**
- **Filters are declarative and cached**
- **Failures degrade gracefully** (empty plots, warnings, not crashes)

---

## Features

- Interactive embedding views (2D / 3D)
- Gene expression visualisation
- Differential expression (volcano plots)
- Dataset-aware filtering (cluster, condition, sample, cell type, genes)
- Deterministic caching of filtered subsets
- Session metadata import/export
- Structured logging (JSON in production)

---

## Repository Structure

single-cell-browser/
├── sc_browser/
│ ├── config/ # Config loading & validation
│ ├── core/ # Dataset, filtering, views
│ ├── metadata_io/ # Session & figure metadata
│ ├── services/ # Export / persistence logic
│ ├── ui/ # Dash layout & callbacks
│ ├── validation/ # Config & session validation
│ └── logging_config.py
├── tests/
├── Dockerfile
├── requirements.txt
└── README.md

---

## Configuration

### Global Config

The application is driven by a **global config directory**, expected to contain:

config/
├── global.json
└── datasets/
├── dataset_1.json
└── dataset_2.json

`global.json` example:

```json
{
  "app_name": "Single-Cell Browser",
  "enable_dataset_management": true
}
Dataset Config
Each dataset is defined by a JSON file:
{
  "name": "Example Dataset",
  "h5ad_path": "data/example.h5ad",
  "embedding_key": "X_umap",
  "obs_columns": {
    "cluster": "leiden",
    "condition": "condition",
    "sample": "sample_id",
    "cell_type": "cell_type"
  }
}
```

Notes:
Only name, h5ad_path, and embedding_key are required
All obs_columns are optional
Missing metadata does not crash the app — affected views degrade gracefully
Dataflow & Robustness
Dataset Loading
.h5ad files are loaded once at startup
Duplicate obs/var names are automatically made unique (with warnings)
Missing files or invalid configs are logged and skipped
At least one valid dataset must load, or startup fails
Filtering
Filters are represented by a FilterState
All subsetting goes through Dataset.subset_for_state
Filtered datasets are cached by filter signature
Empty subsets are valid and handled safely
Views
Views consume filtered datasets, never raw AnnData
Each view declares what filters it requires
Missing metadata results in empty or informational figures, not exceptions

---

## Running Locally (Development)

### Requirements

- Python 3.11+
- pip or virtualenv
- python -m venv .venv 
- source .venv/bin/activate
- pip install -r requirements.txt

- Run the app (example entrypoint):
- export SC_BROWSER_CONFIG_ROOT=./config
- python app.py
- The app will be available at:
- http://localhost:8050

## Running with Docker (Recommended)

### Build

- docker build -t single-cell-browser .
 - Run
 - docker run \
  -p 8050:8050 \
  -e SC_BROWSER_CONFIG_ROOT=/config \
  -v $(pwd)/config:/config \
  -v $(pwd)/data:/data \
  single-cell-browser

 - Gunicorn is used as the production WSGI server.
Environment Variables
Variable	Purpose
SC_BROWSER_CONFIG_ROOT	Path to config directory
SC_BROWSER_DATA_ROOT	Optional data root override
LOG_LEVEL	Logging level (INFO, DEBUG, etc.)
DEBUG	Enable debug logging
PORT	Server port (default: 8050)


## Logging
Structured JSON logging in production
Human-readable logs in development
All dataset loading, filtering, and view execution is logged
Errors are surfaced in the UI without crashing the server

## Testing

### Run tests with:

- pytest

Tests cover:
- Dataset loading & validation
- Subsetting and caching
- Metadata round-tripping
- Failure cases (missing metadata, empty subsets)


