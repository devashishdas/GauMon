# GauMon

**A desktop monitor for Gaussian electronic-structure calculations**  
**Version:** 7.1.1 · **Platform:** Linux · **Interface:** PyQt5 · **License:** MIT

GauMon monitors Gaussian output, visualizes electronic-structure convergence and calculation progress, and provides input validation and evidence-linked diagnostics. It uses a native desktop window rather than a browser or web server. The monitoring and plotting interface does not require a GPU.

> **Important:** GauMon is an independent monitoring utility. It does **not** include, install, license, or replace Gaussian. Synthetic test cases are not quantum-chemical calculations.

## Quick installation (recommended)

Requirements: Linux with a working graphical desktop, Git, and either Conda or Python 3.10+ with `venv` and `pip`.

```bash
git clone https://github.com/devashishdas/GauMon.git
cd GauMon
bash install.sh
```

The installer lets you choose **Conda** or **Python venv**. It installs a private application copy under `~/.local/share/gaumon/app`, validates dependencies and runs tests before finalizing the installation. It optionally creates a desktop-menu entry and, with permission, adds `~/.local/bin` to `~/.bashrc` after backing up the original file. **No `sudo` is required.**

For a noninteractive Conda installation:

```bash
bash install.sh --mode conda --yes --bashrc --desktop
```

For a noninteractive Python venv installation:

```bash
bash install.sh --mode venv --yes --bashrc --desktop
```

To leave `~/.bashrc` untouched, use `--no-bashrc` instead of `--bashrc` and add `~/.local/bin` to your current shell's `PATH` manually.

After installation, open a new terminal or run:

```bash
export PATH="$HOME/.local/bin:$PATH"
hash -r

gaumon --version
gaumon --env-check
gaumon --test
gaumon
```

`gaumon`, `Gaumon`, and `GauMon` are equivalent launcher names and work from any directory **for the installing Linux user**. The GUI starts as a normal maximized window with standard desktop decorations. It does not start a synthetic calculation automatically. The old `--demo-ui` option is removed.

**Conda environment detail:** The wizard uses the dedicated prefix `~/.local/share/gaumon/conda-env`, not the old named environment `gaumon-cpu`. The Python venv alternative uses `~/.local/share/gaumon/venv`. The installed launchers select their own Python environment; users do not need to activate Conda.

See [INSTALL_LINUX.md](INSTALL_LINUX.md) for installation options, safety behavior, troubleshooting, and desktop setup.

## Manual source checkout (alternative)

If you prefer to run GauMon directly from the source tree instead of installing the launcher:

```bash
conda env create -f environment.yml
conda activate gaumon-cpu
python run.py --env-check
python run.py --test
python run.py --gui
```

For an existing `gaumon-cpu` environment:

```bash
conda env update -n gaumon-cpu -f environment.yml --prune
```

This manual path is distinct from the recommended installer. Alternatively, create an isolated Python 3.10+ venv and install `requirements.txt` inside it.

## Main features

| Workspace | Purpose |
|---|---|
| Overview | Current calculation status, key metrics, and four summary plots |
| Run calculation | Explicitly launch a configured Gaussian job or attach to an existing output log |
| Timeline & ETA | Optimization-step and link history; conservative forecasts based on available timing evidence |
| Analysis | SCF, optimization, frequency, and geometry-related summaries |
| Geometry | Lightweight two-dimensional projection of parsed coordinates; **not** a 3D molecular viewer |
| Live log | Follow Gaussian output and inspect evidence-linked diagnostic messages |
| Test cases | Seven Gaussian input templates and synthetic, license-free test fixtures |
| Time graphs | Reported link elapsed/CPU times and call history |
| Input checker/editor | Validate and edit `.gjf`/`.com` input, including route and input-format warnings |
| Advanced time estimator | Conservative step-duration extrapolation without invented historical timestamps |
| Warnings & errors | Source-linked errors, line numbers, and diagnostic guidance |
| Session jobs | Locally saved job references; reattachment is explicit, not automatic |

Additional on-demand scientific workspaces include **Deep SCF**, **Deep QCNR/LinEq1**, and **Deep Frequencies**. SCF iterations, optimization steps, quadratic-convergence linear-equation iterations, and vibrational modes are distinct quantities and are displayed separately. See [DEEP_SCF_GUIDE.md](DEEP_SCF_GUIDE.md) and [DEEP_QCNR_AND_FREQ.md](DEEP_QCNR_AND_FREQ.md).

### Timing and process safety

ETA values are provisional estimates, not guaranteed completion times. GauMon does not invent durations for historical optimization steps lacking usable timestamps. A stale log does not prove a process is hung. An attached external calculation is read-only: **Stop/Force Stop applies only to a process group started and owned by GauMon**, never to an arbitrary external Gaussian process. Closing the monitor does not automatically terminate a running Gaussian job.

## Command-line examples

Commands below work after wizard installation; when running from a source checkout, replace `gaumon` with `python run.py` in an activated environment.

```bash
gaumon                              # native maximized GUI
gaumon --env-check                  # report runtime dependencies
gaumon --test                       # automated tests
gaumon --list-cases                 # available test templates
gaumon --attach /path/to/job.log    # open an existing output log
gaumon --inspect /path/to/job.log   # headless log inspection
gaumon --diagnose /path/to/job.log  # error diagnostics
gaumon --validate /path/to/job.gjf  # Gaussian input validation
gaumon --export-cases ./test_inputs # export sample inputs
```

Only use explicit real Gaussian execution after configuring your licensed Gaussian executable and its environment. GauMon does not configure Gaussian, its scratch directories, or licensing.

## Included scientific test cases

| Key | Molecule / calculation | Notes |
|---|---|---|
| `sp` | H₂ single-point energy | Basic SCF and termination parsing |
| `optfreq` | H₂O optimization + frequencies | Optimization and vibrational records |
| `freq` | H₂O frequencies | Frequency-only workflow |
| `ts` | NH₃ inversion transition-state example | Imaginary-mode interpretation requires checking the displacement |
| `trimer` | Water trimer optimization + frequencies | Illustrative hydrogen-bonded system |
| `stretched` | Stretched H₂ unrestricted SCF | SCF stress case; examine spin contamination |
| `fe2s2` | Fe–S toy single-point cluster | **Not** a validated broken-symmetry or spin-coupling model |

All templates are educational examples rather than benchmark-certified results. Mock Gaussian logs and SCF stress cases exercise parsing without running Gaussian or producing scientifically meaningful energies.

## Development and validation

```bash
python run.py --test
python run.py --env-check
python -m unittest discover -s tests -v
```

The November 2026 local installation log reported **113 passing automated tests**. This is a report of that local run, **not** a guarantee of passing tests on every platform. Graphical rendering, desktop integration, accessibility, and actual Gaussian execution require independent validation on the target workstation. See [RELEASE_VALIDATION.md](RELEASE_VALIDATION.md) and [PUBLICATION_READINESS.md](PUBLICATION_READINESS.md).

For graphical checks on a supported desktop, run `python run.py --gui`. Headless Qt tests are not a substitute for inspecting the real window.

## Known limitations

- Linux is the documented installation target; Windows and native Wayland behavior are not established by the provided tests.
- The geometry workspace is a 2D projection, not a rotatable 3D molecular renderer.
- The parser does not cover every Gaussian link, ONIOM, or Link1 workflow.
- Saved sessions retain explicit file references, but do not automatically resume or acquire ownership of external processes.
- Historical logs may lack sufficient timing information for an ETA.
- Scientific interpretation requires inspection of original Gaussian output; mock data are synthetic.

## Repository documentation

- [INSTALL_LINUX.md](INSTALL_LINUX.md) — Conda/venv installation and launcher behavior
- [INSTALL_QUICK.md](INSTALL_QUICK.md) — manual and quick-start notes
- [DEEP_SCF_GUIDE.md](DEEP_SCF_GUIDE.md) — SCF parsing and plots
- [DEEP_QCNR_AND_FREQ.md](DEEP_QCNR_AND_FREQ.md) — quadratic-convergence and vibrational analysis
- [RELEASE_VALIDATION.md](RELEASE_VALIDATION.md) — regression evidence and remaining checks
- [CONTRIBUTING.md](CONTRIBUTING.md) — development and contributions
- [CITATION.cff](CITATION.cff) — software citation metadata

## Citation and license

GauMon is maintained at [github.com/devashishdas/GauMon](https://github.com/devashishdas/GauMon). Cite the version used in your work; citation metadata are provided in [CITATION.cff](CITATION.cff). See [LICENSE](LICENSE) for the project's license. GauMon does not distribute Gaussian software or proprietary Gaussian output.
