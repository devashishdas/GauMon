#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
for case in sp optfreq freq ts; do
  printf '\n=== Synthetic test: %s ===\n' "$case"
  ./run.sh --mock-case "$case" --timeout 30
done
printf '\nAll four Gaussian-style simulations completed. This does NOT run real quantum chemistry.\n'
