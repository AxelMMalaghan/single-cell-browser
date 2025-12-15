# syntax=docker/dockerfile:1.6
FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PORT=8050

WORKDIR /app

# Runtime deps
RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    libgl1 \
    libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*

# Create a non-root user/group
RUN groupadd --system app && useradd --system --create-home --gid app app

# Install Python deps as root (fine), then drop privileges
COPY requirements.txt /app/requirements.txt
RUN --mount=type=cache,target=/root/.cache/pip \
    pip install --upgrade pip \
    && pip install -r /app/requirements.txt

# Copy code
COPY . /app

# Make sure the non-root user can read the app and write to /app (or better: a subdir)
RUN chown -R app:app /app

USER app

EXPOSE 8050
CMD ["sh", "-c", "gunicorn -b 0.0.0.0:${PORT} app:server --access-logfile - --error-logfile -"]