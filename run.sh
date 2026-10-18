#!/usr/bin/env bash
# Preserve the calling directory for Gaussian jobs and log pickers.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "${CONDA_DEFAULT_ENV:-}" == "gaumon-cpu" ]]; then
  exec python "$ROOT/run.py" "$@"
fi
if command -v conda >/dev/null 2>&1; then
  exec conda run --no-capture-output -n gaumon-cpu python "$ROOT/run.py" "$@"
fi
printf 'ERROR: Conda is not available and gaumon-cpu is not activated.\n' >&2
printf 'Repair Qt/SIP: conda env update -n gaumon-cpu -f "%s/environment.yml" --prune\n' "$ROOT" >&2
exit 2
