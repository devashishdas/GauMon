# Publish GauMon to GitHub

Repository: https://github.com/devashishdas/GauMon

The prepared files are source-only; original user calculations, Gaussian proprietary binaries, checkpoints, outputs, credentials, and local session databases are excluded. Synthetic regression fixtures under `examples/` are intentionally included.

**From the extracted source folder:**

```bash
bash publish_to_github.sh
```

Requires Git and write authentication to the destination repository. Git may ask you to set an author identity with `git config --local user.name 'Your Name'` and `git config --local user.email 'you@example.com'` (replace the examples with your actual values). HTTPS push requires a GitHub-supported credential helper or token. Alternatively switch the remote to `git@github.com:devashishdas/GauMon.git` if SSH auth is configured.

The script initializes `main` if needed, adds `origin` only when missing, refuses a different remote, makes a commit, and performs a **non-force push**. If the remote contains commits, fetch and merge/rebase deliberately; do not use `--force` without inspecting changes.

After push:

```bash
git status
git remote -v
```

Review https://github.com/devashishdas/GauMon/actions to check GitHub Actions tests. Those online CI checks have not been verified in this environment.
