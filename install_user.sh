#!/usr/bin/env bash
# Compatibility alias: route existing install_user.sh users to new wizard.
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec bash "$ROOT/install.sh" "$@"
