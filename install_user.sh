#!/usr/bin/env bash
# Install a per-user GauMon launcher/desktop entry without root or a browser.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ "${CONDA_DEFAULT_ENV:-}" == gaumon-cpu ]]; then
  PYTHON_BIN="$(python -c 'import sys; print(sys.executable)')"
else
  if command -v conda >/dev/null 2>&1; then
    CONDA_BIN="$(command -v conda)"
  elif [[ -x "$HOME/miniconda3/bin/conda" ]]; then
    CONDA_BIN="$HOME/miniconda3/bin/conda"
  elif [[ -x "$HOME/anaconda3/bin/conda" ]]; then
    CONDA_BIN="$HOME/anaconda3/bin/conda"
  elif [[ -x /usr/bin/conda ]]; then
    CONDA_BIN=/usr/bin/conda
  else
    echo 'ERROR: Activate gaumon-cpu before installing, or make conda available.' >&2
    exit 2
  fi
  PYTHON_BIN="$("$CONDA_BIN" run --no-capture-output -n gaumon-cpu python -c 'import sys; print(sys.executable)')"
fi
if [[ ! -x "$PYTHON_BIN" ]]; then
  echo 'ERROR: The selected gaumon-cpu Python is not executable.' >&2
  exit 2
fi
"$PYTHON_BIN" "$ROOT/run.py" --env-check || {
  echo 'ERROR: Install dependencies from environment.yml first.' >&2
  exit 2
}
mkdir -p "$HOME/.local/bin" "$HOME/.local/share/applications"
LAUNCHER="$HOME/.local/bin/gaumon"
if [[ -e "$LAUNCHER" ]] && ! grep -Eq 'GAUMON_USER_LAUNCHER|GAUMON_DESKTOP_INSTALLED' "$LAUNCHER" 2>/dev/null; then
  echo "ERROR: $LAUNCHER exists and was not created by this installer; refusing overwrite." >&2
  exit 3
fi
# Validate the desktop path BEFORE touching an existing command wrapper.
DESKTOP="$HOME/.local/share/applications/gaumon.desktop"
if [[ -e "$DESKTOP" ]] && ! grep -q 'X-GauMon-Installed=true' "$DESKTOP" 2>/dev/null; then
  echo "ERROR: $DESKTOP exists and was not created by GauMon; refusing overwrite." >&2
  exit 3
fi
{
  printf '#!/usr/bin/env bash\n# GAUMON_USER_LAUNCHER\nset -euo pipefail\n'
  printf 'exec %q %q "$@"\n' "$PYTHON_BIN" "$ROOT/run.py"
} > "$LAUNCHER"
chmod 755 "$LAUNCHER"
# XDG desktop Exec argument quoting; no proprietary assets or new logo.
printf '[Desktop Entry]\nType=Application\nName=GauMon\nComment=Gaussian calculation monitoring\nExec="%s" --gui\nTerminal=false\nCategories=Science;Education;\nStartupNotify=true\nX-GauMon-Installed=true\n' "$LAUNCHER" > "$DESKTOP"
chmod 644 "$DESKTOP"
echo "Installed command: $LAUNCHER"
echo "Installed desktop entry: $DESKTOP"
echo 'Run from anywhere: gaumon --gui'
echo 'Or: gaumon --attach /absolute/path/job.log'
echo 'Ensure $HOME/.local/bin is present in your PATH (log out/in if needed).'
