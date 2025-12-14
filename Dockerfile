# syntax=docker/dockerfile:1.6
FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PORT=8050

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    curl \
    libhdf5-dev \
    libgl1 \
    libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /app/requirements.txt

# BuildKit cache speeds rebuilds massively
RUN --mount=type=cache,target=/root/.cache/pip \
    pip install --upgrade pip \
    && pip install -r /app/requirements.txt

COPY . /app

EXPOSE 8050

# Prefer gunicorn for deployment (requires gunicorn in requirements.txt)
# app:server assumes app.py defines `server = app.server`
CMD ["sh", "-c", "gunicorn -b 0.0.0.0:${PORT} app:server"]
