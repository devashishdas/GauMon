#!/usr/bin/env bash
# Publish the extracted GauMon source tree to the user's GitHub repository.
# No force pushes, no credentials in this script, no commit if nothing staged.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
URL='https://github.com/devashishdas/GauMon.git'
if [[ -z "$(command -v git || true)" ]]; then
    printf 'ERROR: git not installed\n' >&2; exit 2
fi
if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    git init -b main
fi
branch="$(git branch --show-current)"
if [[ "$branch" != 'main' ]]; then
    printf 'ERROR: current branch is %s; expected main. Exiting without changing history.\n' "$branch" >&2
    exit 2
fi
if git remote get-url origin >/dev/null 2>&1; then
    existing="$(git remote get-url origin)"
    if [[ "$existing" != "$URL" && "$existing" != 'git@github.com:devashishdas/GauMon.git' ]]; then
        printf 'ERROR: origin currently points to %s; refusing to change it.\n' "$existing" >&2
        exit 2
    fi
else
    git remote add origin "$URL"
fi
printf 'Repository: %s\nBranch: %s\n' "$URL" "$branch"
printf 'Checking files before commit...\n'
git add -A
printf '\nFiles staged:\n'
git diff --cached --stat
if git diff --cached --quiet; then
    printf 'No new changes to commit.\n'
else
    printf '\nCreating source commit...\n'
    git commit -m 'Publish GauMon Desktop v7.1.1 source and tests'
fi
printf '\nUploading with non-forced push to GitHub...\n'
# If main has changed remotely, this push intentionally fails without overwriting it.
git push -u origin main
printf 'Published: https://github.com/devashishdas/GauMon\n'
