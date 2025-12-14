# syntax=docker/dockerfile:1.6
FROM python:3.11-slim

# Prevent Python from writing .pyc files and enable unbuffered logs
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PORT=8050

WORKDIR /app

# --- System deps ---
# Notes:
# - libhdf5-dev is often needed for AnnData/h5py builds
# - libgl1 + libglib2.0-0 helps with matplotlib/umap/scanpy native bits and some rendering deps
# - build-essential needed for any wheels that fall back to source builds
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    git \
    curl \
    libhdf5-dev \
    libgl1 \
    libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*

# --- Python deps ---
# Copy requirements first to leverage Docker layer caching
COPY requirements.txt /app/requirements.txt
RUN pip install --upgrade pip && pip install -r /app/requirements.txt

# If you rely on Plotly static image export, you usually want kaleido installed.
# (If it's already in requirements.txt, this is harmless.)
RUN pip install --no-cache-dir kaleido

# --- App code ---
COPY . /app

# Dash needs to bind to 0.0.0.0 inside Docker
EXPOSE 8050

# Default: run app.py. Change if your entrypoint differs.
CMD ["python", "app.py"]
