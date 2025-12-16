# Single-Cell Browser

![Status](https://img.shields.io/badge/Status-Beta-orange)
![Python](https://img.shields.io/badge/Python-3.9%2B-blue)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED)
![Framework](https://img.shields.io/badge/Framework-Dash%20%7C%20Scanpy-success)

A high-performance, interactive web application for visualizing and exploring single-cell RNA sequencing (scRNA-seq) data. Built on **Dash** and **Scanpy**, it allows researchers to investigate gene expression, differential expression, and cluster composition without writing code.

---

## 📋 Table of Contents

- [Key Features](#-key-features)
    - [Visualization](#visualization)
    - [Data Management](#data-management)
- [Architecture & Data Flow](#-architecture--data-flow)
    - [Directory Structure](#directory-structure)
    - [Data Flow Lifecycle](#data-flow-lifecycle)
- [Prerequisites](#-prerequisites)
- [Installation & Local Development](#--installation--local-development)
    - [Environment Variables](#environment-variables)
- [Configuration & Data Loading](#-configuration--data-loading)
    - [Adding a New Dataset](#adding-a-new-dataset)
- [Docker Deployment (Production)](#-docker-deployment-production)
- [Security](#-security)
- [Testing](#-testing)
- [Code Quality Enforcements](#-code-quality-enforcements)
- [Troubleshooting](#-troubleshooting)
- [License](#-license)
---

## Key Features

### Visualization
* **Dimensionality Reduction:** Interactive UMAP/t-SNE plots colored by cluster, condition, or gene expression.
* **Gene Expression:** Visualise expression levels of specific genes across the entire dataset.
* **Dot Plots:** Compare the percentage of cells expressing genes and their average expression across clusters.
* **Heat Maps:** High-density expression visualization for gene lists.
* **Volcano Plots:** Real-time **Differential Expression (DE)** analysis to find markers between selected groups.

### Data Management
* **"Backed" Mode:** Efficiently handles large datasets using Scanpy's memory-mapping (`backed='r'`), keeping RAM usage low.
* **Dynamic Filtering:** Filter cells by cluster, sample, metadata, or gene expression thresholds in real-time.
* **Session Management:** Save and load analysis sessions (selected genes, filters, view settings).
* **Dataset Import:** Upload new `.h5ad` files directly through the UI (configurable).

---

## Architecture & Data Flow

The application follows a **Hexagonal Architecture** (Ports & Adapters) to separate core logic from the UI.

### Directory Structure
* `sc_browser/core`: Domain logic (Dataset abstraction, Filtering logic).
* `sc_browser/ui`: Dash layouts, components, and callbacks.
* `sc_browser/services`: Data access, file storage, and session persistence.
* `sc_browser/views`: Visualization logic (factory pattern for different plot types).
* `config/`: JSON configurations for datasets.

### Data Flow Lifecycle

1.  **Initialization**:
    * The app reads `config/global.json` and scans `config/datasets/`.
    * Datasets are loaded in **Read-Only Backed Mode** (`.X` remains on disk).
    
2.  **User Interaction** (e.g., Selecting a cluster filter):
    * **Input:** User clicks a checkbox in the UI.
    * **Callback:** Dash triggers a callback in `sc_browser/ui/callbacks/`.
    * **State Update:** The `FilterState` object captures the constraint.
    
3.  **Data Processing**:
    * The `Dataset` core service applies the `FilterState` to the AnnData object.
    * *Optimization:* Instead of copying data, it creates lightweight "Views" of the AnnData object using boolean masks.

4.  **Rendering**:
    * The filtered data is passed to a specific `View` (e.g., `VolcanoPlotView`).
    * The View calculates statistics (if needed) and returns a Plotly Figure.
    * **Output:** The graph updates in the browser.

---

## Prerequisites

* **Python:** 3.9+
* **Data:** Single-cell data must be in `.h5ad` (AnnData) format.
* **Memory:** Sufficient RAM to hold the *index* of your largest dataset (actual data is read from disk).

---

##  Installation & Local Development

### 1. Clone the Repository
```bash
git clone https://github.com/AxelMMalaghan/single-cell-browser.git
cd single-cell-browser
```

### 2. Create Virtual Environment

```
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

### 3. Install Dependencies

```
pip install -r requirements.txt
pip install -r requirements-dev.txt  # For testing/linting
```

### 4. Run the App

```
# Debug mode can be toggled via environment variables
export DEBUG=1 
python app.py
```

Access the app at http://localhost:8050


### Environment Variables

| Variable | Default | Description |
| :--- | :--- | :--- |
| `DEBUG` | `False` | Set to `1` or `True` to enable Dash debug mode (Do not use in production). |
| `PORT` | `8050` | The port the application binds to. |
| `HOST` | `0.0.0.0` | The network interface to bind to. |

---

## Configuration & Data Loading

The application does not assume specific column names in your `.h5ad` files. You map them using JSON configurations.


### Adding a New Dataset

1.  **Upload Data:** Place your `.h5ad` file in the `data/` directory (or the mounted Docker volume).
2.  **Create Config:** Create a JSON config file in `config/datasets/` (e.g., `lung_cancer.json`).

**Configuration Schema:**

```json
{
  "name": "Lung Cancer Study",
  "path": "data/lung_cancer_v1.h5ad",
  "raw": {
    "group": "Oncology",
    "embedding_key": "X_umap" 
  },
  "obs_columns": {
    "cell_id": "index",        
    "cluster": "leiden_0.5",   
    "condition": "treatment",
    "sample": "patient_id",
    "cell_type": "predicted_cell_type"
  }
}
```

## Docker Deployment (Production)

1. Build and Run

```
docker-compose up --build -d
```

2. Volume Mounting

To ensure your data and configurations persist across container restarts, the `docker-compose.yml` mounts the following:

- `./data:/app/data`: Store your `.h5ad` files here
- `./config:/app/config`: Store your JSON configurations here

3. Resource Limits

For production, you should set memory limits in `docker-compose.yml` to prevent upload handlers from crashing the host:

```YAML
services:
  app:
    deploy:
      resources:
        limits:
          memory: 4G 
```

## Security

This application does not have built-in authentication, as it is designed to sit behind a reverse-proxy and only for in-house use.


## Testing

This app uses `pytest` for unit and integration testing

```Bash
# Run all tests
pytest

# Run with coverage report
pytest --cov=sc_browser
```

## Code Quality Enforcements

Code quality is enforced using `ruff`, `black` and `pyright`.

```Bash
# Linting
ruff check . 

# Apply Formatting
black . 

# Type Checking
pyright
```

## Troubleshooting

**1. "Memory" or Container Crash during Upload**
* **Cause:** The uploaded file exceeds the container's available RAM during the Base64 decoding process.
* **Solution:** Increase the Docker memory limit in `docker-compose.yml` or manually place large `.h5ad` files into the `./data` folder instead of using the web uploader.

**2. "KeyError: 'X_umap'"**
* **Cause:** The dataset config specifies `embedding_key: "X_umap"`, but that key does not exist in `adata.obsm`.
* **Solution:** Check your `.h5ad` file keys or update the JSON config to match your data (e.g., `X_tsne`).

**3. Plots are slow to render**
* **Cause:** The dataset might be extremely large, and the `backed='r'` mode is limited by disk I/O speed.
* **Solution:** Ensure the host machine has a fast SSD.

## License 

TBA

