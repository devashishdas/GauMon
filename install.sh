#!/usr/bin/env bash
# Run the standard-library Linux installer with a supported Python (no sudo).
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "${GAUMON_INSTALL_PYTHON:-}" ]]; then
  PYTHON_CANDIDATES=("$GAUMON_INSTALL_PYTHON")
else
  PYTHON_CANDIDATES=(python3.11 python3.12 python3.10 python3 python3.13)
fi
for candidate in "${PYTHON_CANDIDATES[@]}"; do
  if command -v "$candidate" >/dev/null 2>&1 && \
     "$candidate" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)' >/dev/null 2>&1; then
    exec "$candidate" "$ROOT/install.py" "$@"
  fi
done
printf 'ERROR: GauMon installer requires Python 3.10+ on PATH.\n' >&2
printf 'Use GAUMON_INSTALL_PYTHON=/path/to/python3.11 bash install.sh ...\n' >&2
exit 2
