#!/usr/bin/env bash
set -euo pipefail

# Run from the repository so relative resource and upload paths stay valid.
cd "$(dirname "$0")/.."
export HCR_CLEANUP_ENABLED="${HCR_CLEANUP_ENABLED:-1}"
exec "${HCR_PYTHON:-python3}" -m flask --app app --no-debug run \
  --host "${HCR_HOST:-127.0.0.1}" --port "${HCR_PORT:-6789}" \
  --no-debugger --no-reload
