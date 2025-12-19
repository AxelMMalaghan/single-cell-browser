#!/usr/bin/env bash
set -euo pipefail

APP_DIR="{{ app_dir | default('/opt/app') }}"
DASHY_PORT="{{ service_port | default('8000') }}"
APP_PATH="${APP_DIR}/{{entrypoint}}"

exec $APP_PATH