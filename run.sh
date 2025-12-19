#!/bin/bash

# 1. Export the platform-assigned port so the app can see it
export PORT={{service_port}}

echo "Starting application on port: $PORT"

# 2. Start the server using Gunicorn (matching the Dockerfile behavior)
# This binds to the dynamic port and identifies 'server' from 'app.py'
exec gunicorn -b 0.0.0.0:${PORT} app:server --access-logfile - --error-logfile -