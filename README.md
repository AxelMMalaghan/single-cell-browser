# Single-Cell Browser  
*A modular, extensible Dash application for exploring single-cell RNA-seq datasets.*

---

## Overview

The Single-Cell Browser is an interactive, configurable web application designed for scientists who want to explore single-cell datasets **without needing to write code**. The system provides:

- A clean, intuitive user interface  
- Consistent dataset handling through a unified adapter  
- Modular visualisation “views” that plug into the UI automatically  
- Fast filtering and on-demand computation  
- Support for multiple datasets with minimal setup  

The backend is built using:

- **AnnData / Scanpy** for data access  
- **Dash + Plotly** for interactive visualisation  
- A **View Registry** pattern that cleanly separates UI from analysis logic  
- A **Dataset abstraction** that standardises `.obs` metadata across datasets  

---

## Features

### Core
- Load multiple datasets from a JSON config file  
- Standardised mapping of raw `.obs` columns → semantic names (`cluster`, `condition`, `sample`)  
- Interactive filtering by:
  - Clusters  
  - Conditions  
  - Samples  
  - Genes  
- Fully in-memory performance for responsive exploration  

### Visualisation Views
Each view is a self-contained analysis/plotting module:

- **ClusterView** – UMAP/embedding coloured by cluster or metadata  
- **ExpressionView** – Gene expression visualisation  
- **FeatureCountView** – Total counts vs detected features  
- **DotplotView** – Marker overview across cell groups  
- **HeatmapView** – Expression heatmaps  
- **DEView** – Differential expression testing between groups  

Views automatically appear as tabs in the UI.

### Architecture
- Zero hard-coded UI logic  
- All views register themselves with the `ViewRegistry`  
- Dash callbacks delegate work to views via a unified interface  
- Easily extensible — add new views without touching UI code  

---

## Repository Structure


## Data Flow

The app is pointed at a config root which is a directory that contains `global.json` as well as `datasets/` which is a directory that contains the loaded datasets.
`DatasetConfig` then decodes the path where the `.h5ad` files live & obs their respective obs columns mapping.
`from_config` then turns DatasetConfig into the internal `Dataset` object which is responsible for storing the data, normalising the obs_columnn and also precomputes the series for faster filtering.
It also initialises the subset cache:
    `{ (clusters, conditions, samples, cell_types) -> Dataset }`
And the expression cache:
    `{ (gene tuple) -> DataFrame}`

From there, the Dataset objects are then used by Dash to create the views.

### Per-view Data Flow

`compute_data()` creates a pandas dataframe

---

## Installation

### 1. Create and activate a virtual environment

```bash
python3 -m venv .venv
source .venv/bin/activate


pip install -r requirements.txt


pip install pytest black isort


python app.py