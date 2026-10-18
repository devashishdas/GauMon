#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
./run.sh --test 2>&1 | tee latest_tests.log
