#!/usr/bin/env bash
# GauMon launcher for both a source checkout and a managed user installation.
# NEVER changes the caller's working directory: relative Gaussian paths matter.
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Managed installation: run.sh must work even when no Conda shell is activated,
# or when an unrelated Conda environment is active. Staging is marked before
# tests so this very same path is tested during the installation transaction.
if [[ -f "$ROOT/.gaumon-managed-app" ]]; then
  PREFIX="$(cd -- "$ROOT/.." && pwd)"
  MODE="$(sed -n 's/^mode=//p' "$ROOT/.gaumon-managed-app" | head -n 1)"
  case "$MODE" in
    conda) interpreters=("$PREFIX/conda-env/bin/python") ;;
    venv)  interpreters=("$PREFIX/venv/bin/python") ;;
    *)     interpreters=("$PREFIX/conda-env/bin/python" "$PREFIX/venv/bin/python") ;;
  esac
  for interpreter in "${interpreters[@]}"; do
    if [[ -x "$interpreter" ]]; then
      unset PYTHONPATH PYTHONHOME CONDA_PREFIX CONDA_DEFAULT_ENV QT_PLUGIN_PATH QT_QPA_PLATFORM_PLUGIN_PATH
      export PYTHONNOUSERSITE=1
      exec "$interpreter" -s "$ROOT/run.py" "$@"
    fi
  done
  printf 'ERROR: GauMon managed environment missing beneath %s.\n' "$PREFIX" >&2
  exit 2
fi

# Developer/source checkout: honor the explicitly activated environment.
if [[ -n "${VIRTUAL_ENV:-}" ]] || [[ -n "${CONDA_PREFIX:-}" ]]; then
  exec python "$ROOT/run.py" "$@"
fi
if [[ -x "$ROOT/.venv/bin/python" ]]; then
  exec "$ROOT/.venv/bin/python" "$ROOT/run.py" "$@"
fi
if command -v conda >/dev/null 2>&1 && conda env list 2>/dev/null | grep -Eq '^gaumon-cpu[[:space:]]'; then
  exec conda run --no-capture-output -n gaumon-cpu python "$ROOT/run.py" "$@"
fi
exec python3 "$ROOT/run.py" "$@"
