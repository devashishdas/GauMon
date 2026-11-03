# GauMon v5 installer: repaired prefix-environment validation

Run from the directory **containing** the `GauMon/` repository (the patch can
be its sibling); do not move or reclone the existing repository:

```bash
cd /home/leelabgpu6/Software/GAUMON
unzip -o ~/Downloads/GauMon_Easy_Installer_v5.zip
bash gaumon_v5/setup.sh --mode conda --yes --bashrc --desktop
export PATH="$HOME/.local/bin:$PATH"
hash -r
gaumon --version
gaumon --env-check
gaumon --test
gaumon
```

To avoid writing `~/.bashrc`, replace `--bashrc` with `--no-bashrc`.
The installer backs up a regular `~/.bashrc`, rejects symlinks, and adds PATH
only once. No `sudo`, `conda activate`, or system-wide user-account change is
required: this installs a command usable **from any directory by this user**.

---

## Installer location (v4)

The wizard automatically detects either of these layouts:

```text
GAUMON/GauMon/                    # source checkout
GAUMON/gaumon_v4/                 # patch extracted beside checkout
```

or the patch extracted inside the checkout as `GauMon/gaumon_v4/`.
From `GAUMON/`, run `bash gaumon_v4/setup.sh`; no moving directories
or second clone is needed. The wizard does not touch `~/.bashrc` without
approval (or the explicit `--bashrc` flag).

# GauMon — Easy Linux installation (maximized GUI)

**From an existing checkout, run `bash install.sh` and select Conda or Python venv.**
The wizard asks before touching `~/.bashrc`, backs up the existing file, and
adds an idempotent guarded PATH block only with permission. The application
opens as a **normal maximized window** (title bar and taskbar visible), never full-screen.
`--demo-ui` is no longer a supported command; synthetic mock tests remain
available through `--test` and `--mock-case`.

After opening a new terminal, type `gaumon`, `Gaumon`, or `GauMon` from any directory.
This is per-user availability, not a root-installed multi-account executable.

# GauMon v7.1.1 — Linux installation (Conda and non-Conda)

GauMon is a graphical Gaussian log monitor. It does **not** ship or install Gaussian, invoke Gaussian without an explicit user action, or alter existing jobs. Installation is per-user; **do not use sudo**.

## Supported setups

- Linux x86_64, X11 or a working XWayland desktop (including Rocky Linux / XRDP). This installer does not yet claim tested Wayland-native or Windows support.
- Conda path: Miniforge/Miniconda/Anaconda with internet access and the `conda-forge` channel.
- Without Conda: Python **3.10+** with `venv`/`ensurepip`, `pip`, and internet access to PyPI. Python **3.11** is recommended for reproducibility.
- A Linux GUI requires functional Qt/XCB platform libraries. Virtual/headless computers can still use analysis commands.

### New machine quick start

```bash
sudo dnf install -y git python3 python3-pip  # Rocky/Fedora, only if missing
# On Ubuntu/Debian instead: sudo apt update && sudo apt install -y git python3 python3-venv python3-pip

git clone https://github.com/devashishdas/GauMon.git
cd GauMon
bash install.sh
```

The wizard offers **1. Conda** and **2. Python venv** and asks whether to add `~/.local/bin` to `~/.bashrc` and install an application menu shortcut. It installs a **stable source copy** at `~/.local/share/gaumon/app`, rather than depending on the Git checkout location.

### Fully noninteractive Conda

```bash
bash install.sh --mode conda --yes --bashrc --desktop
```

Environment: `~/.local/share/gaumon/conda-env/`. The wizard creates it from `environment.yml` (or updates it in place), then runs `--env-check`, `--list-cases`, and the unit tests. If Conda cannot be found:

```bash
bash install.sh --mode conda --conda-exe "$HOME/miniforge3/bin/conda" --yes --bashrc
```

This uses a dedicated **prefix-based** environment (the old environment **name** `gaumon-cpu` is not required and is not changed).

### Fully noninteractive without Conda

```bash
bash install.sh --mode venv --yes --bashrc --desktop
```

Creates `~/.local/share/gaumon/venv/`; installs dependencies with `venv/bin/python -m pip install -r requirements.txt` and validates with `pip check`. An externally managed OS Python (PEP 668) is **not** written to because all dependencies go into the private venv.

### Use from any working directory

Once the installer has completed and PATH has been refreshed:

```bash
command -v gaumon
gaumon --version
gaumon --env-check
gaumon                 # launches maximized GUI
gaumon --attach /absolute/path/job.log
gaumon --inspect /absolute/path/job.log
gaumon --list-cases
```

If you did **not** choose `--bashrc`, enable the PATH for only your current shell, or opt into a safe persistent edit:

```bash
# Safer alternative: rerun installer with explicit, backed-up .bashrc opt-in:
# bash install.sh --mode conda --yes --bashrc
# For the current shell only, without changing files:
export PATH="$HOME/.local/bin:$PATH"
```

The installer only edits the user's `~/.bashrc` on request. A managed guarded `GauMon PATH` section avoids repeated entries in `PATH` and duplicate additions on repeat installations. Existing `.bashrc` content and permissions are preserved. It backs up an existing shell config as `~/.bashrc.gaumon-backup` before modifying it. The per-user launchers are `gaumon`, `Gaumon`, and `GauMon` in `~/.local/bin`, not `/usr/bin/gaumon`; they are available in every shell where `~/.local/bin` is in `PATH`.

### Dry run and alternate install location

```bash
bash install.sh --mode auto --yes --dry-run
bash install.sh --mode venv --yes --prefix "$HOME/apps/gaumon-install" --no-desktop
```

The installer refuses to overwrite non-GauMon launchers, desktop files, or application directories. Existing managed application content is replaced only **after** dependency checks and tests pass. Re-running the installer updates a matching private environment.

### Runtime GUI preflight

```bash
gaumon --env-check
gaumon --test
gaumon --gui
```

An interactive Qt/XCB smoke test is still needed. On a *new* Linux system where the GUI cannot connect to the display, verify `DISPLAY`, X11 forwarding/XRDP, `libxcb`, and Qt-related system libraries. Useful diagnostics:

```bash
echo "$DISPLAY"
QT_DEBUG_PLUGINS=1 gaumon --gui
```

Do not use `QT_QPA_PLATFORM=offscreen` for normal interactive usage (it is suitable for CI only). If `PyQt5.sip` fails, reinstall the private environment rather than mixing operating-system Python packages or copying `.so` files.

### Development (without installing a launcher)

```bash
# Conda
conda env create -f environment.yml
conda activate gaumon-cpu
python run.py --test
python run.py --gui

# Python venv
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python run.py --test
python run.py --gui
```

The wizard accepts `--skip-tests` for debugging, but checks imports and packaged examples regardless. Since virtual machines and fresh Linux distributions vary, report the OS, Python/Conda version, `gaumon --env-check` output, and any `QT_DEBUG_PLUGINS` error when opening an issue.

## Repair for Rocky Linux: `PermissionError: examples/mock_g16`

If `--test` or `--mock-case sp` reports `PermissionError: [Errno 13]` for
`examples/mock_g16`, the checkout or archive lost its executable bit:

```bash
cd /home/leelabgpu6/Software/GAUMON/GauMon
chmod +x examples/mock_g16
python run.py --test
```

**Installer fixes (2026-11-03):** Installing to a persistent prefix now
copies both `run.sh` and `run_tests.sh`, restores `mock_g16` executable
permission in the staged copy, and validates the staged run script under
its actual private environment before promotion. This addresses the missing
`run.sh` tests and repeated `PermissionError` observed on Rocky Linux.

For publication, also commit the executable bit to Git:

```bash
git add --chmod=+x examples/mock_g16
git ls-files -s examples/mock_g16  # expected 100755 after staging
```

The installer can safely be rerun after a failed staging attempt. It retains
the dedicated environment and will revalidate the corrected application.
Do **not** delete your existing `gaumon-cpu` environment or Gaussian install.
