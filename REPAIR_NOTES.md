## v4: Fix source checkout discovery

`setup.sh` previously assumed its parent was the GauMon checkout.
When the patch and clone were siblings (`GAUMON/GauMon/` and `GAUMON/gaumon_v3/`),
it failed before installation. The setup entrypoint now detects both
sibling and embedded layouts and preserves installer arguments without
touching shell configuration itself.

# GauMon installation repair (2026-11-03)

## Evidence from Rocky Linux / leelabgpu6

- The original checkout passed Qt dependency preflight, but 7/90 tests
  and `--demo-ui` failed with `PermissionError` on `examples/mock_g16`.
- Installer staging (Conda prefix under `~/.local/share/gaumon`) reached
  dependency preflight, then failed 9 of 100 tests: seven `mock_g16`
  execute-permission errors, one missing `run.sh` error, and one failed
  relative-working-directory test due to absent `run.sh`.
- The original activated `gaumon-cpu` leaked `CONDA_PREFIX` and Qt-related
  environment variables into staged preflight. They are now sanitized.

## Permanent code corrections

1. `apply_patch.sh` repairs executable permissions for the checkout's
   `examples/mock_g16` and installs a source/managed-prefix-aware `run.sh`.
2. `install.py` copies `run.sh` and `run_tests.sh` to the managed stage.
3. Every copy explicitly sets execute permission on `examples/mock_g16`,
   `run.sh`, and `run_tests.sh` even if an archive lost POSIX mode bits.
4. Staging validation refuses a missing or non-executable mock launcher.
5. The staged app receives a managed marker *before* testing. This means
   `run.sh` uses the newly created environment even before promotion.
6. The marker stores the selected installer mode so `run.sh` never picks
   an unrelated existing Conda/venv environment during mode changes.
7. Added tests for permission repair, missing `run.sh`, and environment cleanup.

## Apply to an existing checkout

```bash
cd /home/leelabgpu6/Software/GAUMON/GauMon
unzip -o ~/Downloads/GauMon_Maximized_Easy_Installer_v3.zip
bash gaumon_v3/setup.sh --mode conda --yes --bashrc --desktop
python run.py --test
# The one-command setup script also installs/updates the managed copy.
~/.local/bin/gaumon --env-check
~/.local/bin/gaumon --test
~/.local/bin/gaumon --gui
```

No need to delete the existing working `gaumon-cpu` Conda environment.
The wizard uses a separate private prefix, and will safely update it.

## GitHub POSIX executable bit

```bash
git add --chmod=+x examples/mock_g16
git add run.sh run_tests.sh install.py install.sh install_user.sh \
  tests/test_install_wizard.py .github/workflows/ci.yml
```

The executable bit must be committed so *fresh* `git clone` checkouts work,
not just local copies patched with `chmod`. Verify with
`git ls-files -s examples/mock_g16` after staging; expected mode is `100755`.

## Verification status

The new installer has offline regression tests. Actual Linux desktop/Qt
rendering and a fresh machine Conda/venv solve must still be tested on the
user's workstation. This repair does not claim scientific validation of
Gaussian output parsing or readiness for JOSS acceptance.

## v3 behavior
- `--demo-ui` is removed from the public CLI.
- Native Qt window opens maximized; F11/Fullscreen button removed.
- `gaumon`, `Gaumon`, and `GauMon` are available from any directory after PATH initialization.
- `~/.bashrc` edits require an explicit wizard choice or `--bashrc`; existing content is preserved, backed up, and the PATH entry is guarded against duplication.


## v5 (2026-11-03): installer aborts at 111/112 tests

**Observed error:** `test_env_check_reports_expected_name` requires the literal
`gaumon-cpu` in an environment report. This assertion is incompatible with the
new isolated-prefix Conda installation, which uses
`~/.local/share/gaumon/conda-env` and deliberately clears `CONDA_PREFIX` to
avoid environment and Qt contamination. The dependencies themselves passed.

**Fix:** `patch_source.py` now replaces the stale test with assertions that
check `run.py --env-check` exit status and real PyQt5/pyqtgraph imports. It does
not weaken dependency validation, bypass the test suite, or force a fake
`CONDA_PREFIX` environment variable. The patch is idempotent and refuses an
unrecognized source instead of silently changing unrelated tests.

**Install** (archive extracted next to `GauMon/` or inside it):

```bash
cd /home/leelabgpu6/Software/GAUMON
unzip -o ~/Downloads/GauMon_Easy_Installer_v5.zip
bash gaumon_v5/setup.sh --mode conda --yes --bashrc --desktop
export PATH="$HOME/.local/bin:$PATH"
hash -r
gaumon --test
gaumon
```

The installation completes and adds the launchers **only after** `--env-check`,
`--list-cases`, and all unit tests pass. `~/.bashrc` is unchanged unless `--bashrc`
is explicitly selected or the interactive wizard receives permission.
