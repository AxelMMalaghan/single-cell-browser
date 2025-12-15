FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PORT=8050

WORKDIR /app

# Runtime deps (keep lean; add build-essential only if wheels fail)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    libgl1 \
    libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /app/requirements.txt

RUN --mount=type=cache,target=/root/.cache/pip \
    pip install --upgrade pip \
    && pip install -r /app/requirements.txt

COPY . /app

EXPOSE 8050
CMD ["sh", "-c", "gunicorn -b 0.0.0.0:${PORT} app:server"]