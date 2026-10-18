# New in v7.1.1

- **Full screen on launch** (F11 toggle; Esc maximize), with `GAUMON_START_MAXIMIZED=1` XRDP fallback.
- **CWD defaults** for Gaussian new-job directory, file choosers, exports and running logs.
- **Auto input viewer:** starting a `.gjf` job displays that exact source; attaching a `.log` loads a verified matching input in the same folder, or reports when it cannot safely identify the source. Read-only resumed sessions retain explicit input provenance.
- **Keyword guide:** line-numbered input editor now provides explanations for Gaussian `%` / `#` / route keywords and important scientific cautions.
- **Advanced ETA deque crash fixed** and covered by targeted regression tests.

# GauMon — standalone Gaussian monitor

**Standalone Gaussian 16 monitor · native PyQt5 window · no browser · no GPU requirement.**

GauMon Desktop v7.1.1 integrates the original v5.2 deep-analysis engine with the five requested utility pages. It provides a calm, high-contrast-on-white scientific interface inspired by the approved light dashboard concept. It is a **working test release**: monitoring, parsing, safe owned-process control, seven Gaussian input templates, optional Deep SCF, Deep QCNR/LinEq1, and Deep Frequencies workspaces, and automated tests are implemented. Actual **visual inspection of the Qt window on Rocky Linux is still required** before claiming a production-quality release.

## Interface principles

- No logo yet, no blue/black gaming appearance, no pop-up browser, no HTTP server, no CUDA/OpenGL backend.
- Warm neutral off-white, muted forest green, limited coral/orange diagnostic colors.
- Consistent 13 px body type and ≥20 px metric values, fully native OS text rendering.
- Dashboard contains only **four central charts + a small run-details panel**. Setup, table, logs, geometry and tests are separate navigation views rather than crammed into the dashboard.
- **No application-wide scroll area.** The window uses fluid Qt layouts at 1366×768 and 1920×1080 target screen sizes; tables and log text have their own scrollbars because data can be arbitrarily long. The minimum window size is 1160×680. Native display scaling, remote desktop panels and missing fonts can affect appearance.
- Plots use pyqtgraph's **CPU QPainter mode**; the application explicitly sets `useOpenGL=False`.

## Installation: exact Conda environment

**Use the environment named `gaumon-cpu`, not Conda `base`, `boltz2_cuda121`, or a Gaussian-specific environment.** This avoids Qt version conflicts. The `environment.yml` pins Python 3.11 and PyQt 5.15, and explicitly includes `pyqt5-sip`, numpy, pyqtgraph, and psutil from conda-forge. These packages are for monitoring; Gaussian binaries are installed and licensed independently.

```bash
# From the extracted GauMon_Desktop_CPU_v7_1_1 directory
conda env create -f environment.yml
conda activate gaumon-cpu
python run.py --env-check
python run.py --test
python run.py --demo-ui
```

If the environment already exists, update rather than recreate it:

```bash
conda env update -n gaumon-cpu -f environment.yml --prune
conda activate gaumon-cpu
```

For **Rocky Linux 8 / XRDP**, `python run.py --gui` or `./run.sh` works with software Qt rendering. If Qt reports an XCB failure, inspect system `libxcb`/X11 libraries (not GPU drivers) and use `QT_DEBUG_PLUGINS=1 python run.py --gui` to locate the missing shared library. Do not set `QT_QPA_PLATFORM=offscreen` for your normal interactive session; that flag is reserved for CI/headless GUI testing.

`./run.sh` uses `gaumon-cpu` through `conda run` even when the shell is not activated. It refuses to silently use an unrelated Python environment. If your `conda` shell command is unavailable, initialize Conda with your installation's `conda.sh` then retry.

An alternative pip requirements file is provided for people who already have a working Python 3.10+ environment:

```bash
python -m pip install -r requirements.txt
python run.py --gui
```

## Start GauMon

```bash
conda activate gaumon-cpu
python run.py --gui             # main desktop app
python run.py --demo-ui         # opens GUI and starts synthetic opt+freq
python run.py --list-cases      # show shipped Gaussian test inputs
python run.py --env-check       # dependency / environment diagnostics
```

### Fifteen views — seven standard, three deep scientific, five utilities

Optimization-cycle timing, original step/link-event history, and four ETA horizons are available on the dedicated Timeline & ETA page. Historical steps without validated timestamps remain untimed; the application does not assign zero-second durations when it first reads an attached log.


| View | Contents |
|---|---|
| **Overview** | Six live metrics, four plots, current link / SCF cycles, and one compact ETA preview; failures display a small diagnostic banner |
| **Run calculation** | Gaussian executable, `.gjf`/`.com` input, working directory, output filename; Start / Stop gracefully / Force stop; external read-only log attachment |
| **Timeline & ETA** | Four forecast cards for the next +10 / +50 / +100 / +150 optimization cycles, cycle history, and repeated Gaussian link history |
| **Analysis** | Complete SCF history, frequency/optimization summary, and geometric warning thresholds |
| **Geometry** | Simple 2D XY projection of the latest parsed coordinates; it is **not** a full 3D molecular renderer |
| **Live log** | Actual Gaussian output tail and specific, source-credited Gaussian failure diagnosis; log-only scrolling is intentional |
| **Test cases** | Four basic and three advanced scientific input templates; synthetic simulation; export selected/all `.gjf` inputs |

### Five additional utility workspaces

| New page | Purpose |
|---|---|
| **Time graphs** | Measured `Leave Link` elapsed and CPU durations, link totals, call sequence; explains wall vs CPU and why no per-atom time is available |
| **Input checker / editor** | Editable colored `.gjf`/`.com`, line numbers, whitespace/line-ending visualization, problem navigation, and safe Save As |
| **Advanced time estimator** | Observed timing for newly followed optimization steps, conservative horizon estimates; refuses to invent historical durations |
| **Warnings & errors** | Evidence-linked Gaussian warnings with line numbers and nonautomatic diagnostic guidance |
| **Session jobs** | Up to 20 explicitly opened/launched logs; saved local records must be manually resumed read-only after restart |

The original ten pages are preserved. See [MERGE_RELEASE_NOTES.md](MERGE_RELEASE_NOTES.md) for the combined architecture and [INSTALL_QUICK.md](INSTALL_QUICK.md) for setup, including the PyQt5.sip repair procedure.

**Deep workspaces are hidden by default.** On the Overview page choose:

- **Deep SCF diagnostics →**: DIIS/SCF iteration graphs; reset for each new SCF procedure.
- **Deep QCNR / LinEq1 →**: nested linear-equation residuals within quadratic-convergent SCF (L508); reset for each new LinEq1 solve. This is **not vibrational frequency** output.
- **Deep Frequencies →**: actual `Frequencies --` normal modes, reduced masses, force constants, and IR intensities where printed. A single-point run without frequency analysis shows no frequency modes.

Each advanced workspace has a history selector (as applicable), charts, and export controls. Closing it returns to Overview. See [DEEP_QCNR_AND_FREQ.md](DEEP_QCNR_AND_FREQ.md).

**Process state:** `STARTING`, `ATTACHED` (historical log, unverified), `RUNNING`, `STALLED`, `COMPLETED`, `FAILED`, `STOPPED`, or `INCOMPLETE`. `STALLED` means the log has not grown for three minutes; it is not proof of a hung job. The observed timer freezes at a verified terminal state. A historical log attached without a PID shows its observation duration, **not** Gaussian's original runtime.

**Process safety:** GauMon's Stop/Force stop commands operate only on a process group that GauMon itself launched. The program never runs a global `pkill`, never kills an attached external process, and does not terminate Gaussian merely because the GUI closes. Closing the app while a managed process remains active requires confirmation; the job continues afterward.

## Use GauMon from anywhere in the computer

Install a command and desktop application entry for the **current user** without changing the Gaussian installation:

```bash
conda activate gaumon-cpu
./install_user.sh
# Ensure ~/.local/bin is in PATH (open a new shell if needed)
export PATH="$HOME/.local/bin:$PATH"
gaumon --gui
gaumon --attach /absolute/path/job.log
gaumon --inspect /absolute/path/job.log
gaumon --list-jobs
```

The launcher uses the actual Python executable inside `gaumon-cpu` and preserves the caller’s working directory. **Keep this source folder at its installed location** because the launcher references it. The installed desktop shortcut appears as **GauMon** in Linux application menus; it includes no custom logo. Running jobs discovered by PID can be chosen inside **Run calculation** using *Find running Gaussian processes*. Attached jobs are read-only. Without a PID, a previously written log is shown as **ATTACHED** until new output arrives; it is not incorrectly claimed to be running.

## Why the calculation failed: FormBX / L103

The signature `Bend failed for angle ...`, `Tors failed for dihedral ...`, `FormBX had a problem.` with `l103.exe` indicates Gaussian's redundant-internal-coordinate optimization could not proceed. This often occurs near a linear or otherwise ill-defined bond-angle/torsion arrangement; a detailed geometry inspection is needed to establish the actual cause. The monitor extracts the numbered atoms from the log, prioritizes **FormBX** over generic error termination, and shows a compact failure banner and a full explanation under **Live log**. A supplied regression fixture shows the failed angle **19–27–25** and torsions **20–19–27–25**, **21–19–27–25**, **8–25–27–19**, and **26–25–27–19**.

For an unconstrained optimization, consider inspecting the last valid structure and rerunning with `Opt=Cartesian`. This can avoid internal-coordinate problems but may require more iterations. Do not blindly combine it with `Opt=ModRedundant`; review frozen/scan variables and chemistry before changing the route. This guidance is **credited to Zhe Wang, [Gaussian Common Errors and Solutions](https://wongzit.github.io/gaussian-common-errors-and-solutions/)** (accessed 2026-10-18). The application displays the credit and includes a link in the log diagnostics panel. The fixture contains only a representative user-supplied failure excerpt, not a scientifically validated result.

## ETA and cycle timing: correct interpretation

- The parser records every **distinct reported optimizer step** and every **`Enter ... lXXX.exe`** event, even if an existing log was written before you attached it.
- For cycle durations it first attempts to use elapsed differences between timestamp proxies from nearby Gaussian `Leave Link ... at ...` records. If unavailable, it estimates durations from **newly observed live** step changes only. The first step after attachment establishes an observation baseline.
- Multiple historical steps parsed together **do not** become zero-time cycles. If there is insufficient timing data, the ETA shows an em dash instead of fabricated dates.
- For ongoing jobs, the **recent-cycle average** extrapolates +10, +50, +100, and +150 *additional* optimization steps. Times are displayed in the machine's local timezone, and change as new cycle intervals arrive. The final real completion time is not guaranteed: Gaussian may converge early or switch to costly operations such as frequencies.
- Terminal jobs do not display predictions for a calculation that already finished. The elapsed clock measures monitor observation duration; attach-to-existing historical runtime cannot be reconstructed from a log with insufficient time records.
- Link history preserves repeated link entrances (e.g. L103 → L502 → L103), and the provided `gaussian_links.json` supplies link descriptions.

## Gaussian test cases: four fundamentals plus three advanced

All seven templates contain `%mem`, `%nprocshared`, a checkpoint name, route section, **nonempty title card**, blank separators, charge/multiplicity and Cartesian coordinates. None asks for `geom=connectivity` without supplying a connectivity table.

| Key | Input | Molecule | Gaussian job | What to verify on **real** output |
|---|---|---|---|---|
| `sp` | `01_energy_H2.gjf` | H₂ | `#p b3lyp/3-21g sp` | SCF Done and normal termination |
| `optfreq` | `02_opt_freq_H2O.gjf` | H₂O | `#p b3lyp/3-21g opt freq` | Optimization completed, frequency records, normally no imaginary modes |
| `freq` | `03_freq_H2O.gjf` | H₂O | `#p b3lyp/3-21g freq` | Frequencies printed; no optimization required |
| `ts` | `04_ts_NH3_inversion.gjf` | NH₃ (planar) | `#p b3lyp/3-21g opt=(ts,calcfc,noeigentest,maxcycles=35) freq nosymm` | Saddle-point optimization and **one** imaginary frequency along umbrella inversion if it converges |

**Chemistry caveats:** These are short educational examples, not benchmark-certified results. TS convergence and vibrational frequency signs depend on method, geometry and numerical settings. A frequency-only run at an unoptimized structure can produce modes that do not characterize a stationary point. Verify the displacement associated with any imaginary mode before claiming a TS.


### Advanced inputs and SCF stress tests

The **three additional inputs** are illustrative stress tests for monitoring and parsing. They are not guaranteed to converge, nor are the chosen geometries and electronic states chemically validated.

| Key | Input | Purpose | Scientific caveat |
|---|---|---|---|
| `trimer` | `05_advanced_water_trimer_optfreq.gjf` | Larger hydrogen-bonded geometry optimization plus frequency | Alternative minima or rearrangements possible |
| `stretched` | `06_advanced_stretched_H2_UHF.gjf` | Difficult stretched-bond unrestricted SCF with `Guess=Mix`/`XQC` | Examine spin contamination and solution stability |
| `fe2s2` | `07_advanced_Fe2S2_toy_SP.gjf` | Small metal–sulfur electronic-structure stress case | **Not** an experimentally validated Fe–S spin-coupling or broken-symmetry model |

The program also includes four deterministic **synthetic SCF stress logs** (no licensed Gaussian required): `many_steps`, `oscillation`, `incomplete`, and `chunked`. From **Test cases**, choose an SCF stress fixture and click *Inspect in Deep SCF*. You can also use the command line:

```bash
python run.py --stress-case many_steps
python run.py --stress-case oscillation
python run.py --deep-report /path/to/d.log
python run.py --deep-report /path/to/WaterDrops.log
```

### Deep SCF: electronic iterations within each geometry optimization

**SCF cycle ≠ optimization step.** A Gaussian `Cycle N Pass P IDiag I` record reports an electronic SCF iteration inside one SCF procedure, often at a fixed nuclear geometry. The `Step number ...` record refers to a geometry-optimization step, which typically follows an SCF result. One optimization step can sometimes require more than one SCF call.

The **Deep SCF** view is user-activated and contains four plot areas:

1. Electronic energy `E` versus SCF cycle (Hartree).
2. Absolute `Delta-E` versus SCF cycle (Hartree; logarithmic display).
3. DIIS error and `RMSU` versus SCF cycle (diagnostic, logarithmic display).
4. `RMSDP` and `MaxDP` versus SCF cycle (density change, logarithmic display).

The compact iteration table displays cycle, pass, energy, energy difference, DIIS error, RMSU, RMSDP, and MaxDP. Selecting a row reveals further `OVMax`, `NSaved`, gaps, CP values, coefficients (`Coeff`, `Coeff-Com`, `Coeff-En`), weights and DIIS matrix descriptors. Export the selected SCF block as CSV for external analysis.

**Reset policy:** The *Follow newest* selection shows only the latest SCF procedure. On the next `Cycle 1`/cycle-index reset after a completed `SCF Done`, the graphs start with a new block rather than appending to the old curves. Completed SCF blocks are preserved in a bounded history (at most 120 procedures × 400 iterations each), selectable from the dropdown. The view does not reset the whole calculation or delete any Gaussian output.

Partial or unfinished SCF calculations remain identified as *running/incomplete data*. A `SCF Done` marker is required to mark an SCF block complete; an optimization can fail or terminate before reaching that marker. Discarded historical rows are signaled as *truncated*. Zero and missing values are omitted from logarithmic plots, not converted into invented values.

**Real-log checks during development:** `WaterDrops.log` contained 19 complete SCF procedures, with the final `Step number` not yet printed for SCF #19. `d.log` contained one ongoing 19-cycle SCF procedure, without a final `SCF Done`. The raw uploaded Gaussian logs are **not distributed** with the app. Use your own authorized logs locally.

### Safe, license-free testing

```bash
python run.py --test                              # automated regression tests
python run.py --mock-case sp --timeout 30
python run.py --mock-case optfreq --timeout 30
python run.py --mock-case freq --timeout 30
python run.py --mock-case ts --timeout 30
```

A mock process emits **synthetic Gaussian-style text** and **does not perform quantum chemistry**. The `examples/mock_logs/` folder contains fixture logs for the four basic calculation types, and `examples/stress_logs/` contains deep SCF stress fixtures. The synthetic TS has exactly one *fabricated* negative frequency solely to exercise the warning and parsing logic.

### Explicit real Gaussian execution

Do **not** accidentally launch Gaussian from the automated unit tests: real runs are opt-in.

```bash
# Export all seven .gjf inputs, without overwriting existing inputs
python run.py --export-cases ./test_inputs

# Verify input syntax first
python run.py --validate ./test_inputs/02_opt_freq_H2O.gjf

# Only when your Gaussian 16 executable and license are ready:
python run.py --run-real sp --gaussian-exe g16 --output ./real_sp --timeout 300
python run.py --run-real optfreq --gaussian-exe g16 --output ./real_optfreq --timeout 1200
python run.py --run-real freq --gaussian-exe g16 --output ./real_freq --timeout 600
python run.py --run-real ts --gaussian-exe g16 --output ./real_ts --timeout 1200
```

You can use a full executable path, e.g. `--gaussian-exe /home/exampleuser/g16/g16` **if that is the correct executable on your system**. The executable path, module environment (`g16root`, scratch), and license must be configured according to your installation. The `gaumon-cpu` Conda environment does not install Gaussian.

The test-case runner refuses to overwrite an existing input file or nonempty output log. Use a fresh output directory when repeating a calculation.

## Diagnosis: your reported Gaussian error

The supplied Gaussian output contains the sequence:

```
(Enter /home/exampleuser/g16/l101.exe)
Title Card Required
...
End of file in ZSymb.
Error termination via Lnk1e ...
Error: segmentation violation
```

**First explicit parsing error:** `End of file in ZSymb.`, in Gaussian's initial input-reading link `l101.exe`. A missing terminating blank line after the geometry is a strong candidate. `Title Card Required` can also be a literal placeholder title, and by itself is **not** proof of a missing title section. The subsequent segmentation-violation traceback should not be assumed to be the primary cause without first fixing the input.

A corrected **format** looks like this (keep your original atoms, charge and multiplicity; these are not shown here):

```text
%nprocshared=16
%mem=32GB
%chk=Arg_opt_freq.chk
#p opt freq b3lyp/aug-cc-pvtz

Arginine opt and frequency calculation

1 1
<your original Cartesian coordinates here>

```

The angle-bracketed line is explanatory only and **must** be replaced by your original atomic-coordinate lines; it is not runnable as written. In the log you pasted, `geom=connectivity` also appears. **Remove `geom=connectivity` if you have no explicit connectivity block**, or supply a correct block after the coordinates. Do not change charge `1` or multiplicity `1` without reviewing your actual molecular model.

The user-supplied output prints atomic coordinates before `End of file in ZSymb`, which strengthens the case for a missing final separator or incomplete geometry/connectivity input. The original `.gjf` file was not supplied, so this diagnosis remains provisional. For the seven shipped tests, valid input templates are included. You can diagnose a Gaussian output without launching the GUI:

```bash
python run.py --validate /path/to/your.gjf
python run.py --diagnose /path/to/your.log
python run.py --diagnose examples/mock_logs/zsymb_eof_example.log
```

If a minimal validated input still crashes `l101.exe` with segmentation violation, inspect your Gaussian 16 Rev. C.01 installation, its runtime environment, license, and system libraries. An input-format error does **not** establish that GauMon caused a Gaussian executable crash.

## Diagnostics and test expectations

```bash
./run_tests.sh                     # captures headless test log
python run.py --env-check
python run.py --demo-ui             # manual GUI smoke test
```

Expected **synthetic** results (not actual quantum chemistry):

- `sp` → `COMPLETED`, one SCF energy, zero frequency blocks.
- `optfreq` → `COMPLETED`, optimization marker and frequency blocks.
- `freq` → `COMPLETED`, frequency blocks without optimizer-completed marker.
- `ts` → `COMPLETED`, optimizer marker and exactly one negative mock frequency.
- A managed stop → `STOPPED`, frozen observed timer.
- An invalid input lacking a title → validator flags a missing title; an input lacking the final blank line → validator warns about `End of file in ZSymb` risk.

## Project layout

```text
GauMon_Desktop_CPU_v5/
├── run.py                       # CLI + GUI entry point
├── run.sh / run_tests.sh        # use dedicated Conda environment
├── install_user.sh              # per-user command and desktop shortcut
├── environment.yml             # gaumon-cpu (Python 3.11 + PyQt5)
├── requirements.txt            # optional pip alternative
├── gaumon/
│   ├── ui.py                    # native CPU-only dashboard, 7 views
│   ├── engine.py                # incremental output monitoring + safe process groups
│   ├── inputs.py                # Gaussian input templates + validation
│   ├── diagnostics.py           # known log/error diagnoses
│   └── __init__.py
├── examples/
│   ├── inputs/                  # 4 real Gaussian .gjf input templates
│   ├── mock_logs/               # 4 synthetic logs + title-error fixture
│   ├── mock_g16.py / mock_g16   # mock command: no Gaussian installation
├── assets/gaussian_links.json
└── tests/                      # headless regression + optional Qt smoke tests
```

## Known limitations and release gates

1. **GUI rendering must be checked on your Rocky/XRDP workstation.** Qt dependencies are unavailable in this build environment. Python compiles and headless tests run here, but the exact pixel layout has **not** been verified by taking a live Qt screenshot. Run `python run.py --demo-ui` at both 1366×768 and 1920×1080 before publishing the app.
2. No GUI-global scroll area is used. On displays smaller than the 1080×650 minimum, widgets may not fit. Increase desktop resolution or OS workspace size rather than squeezing scientific plots below readable dimensions.
3. Geometry projection is a CPU QPainter 2D view, not rotatable three-dimensional molecular graphics.
4. External PID attachment is user-supplied; file growth without a verified process does not prove Gaussian is currently running.
5. Parsing handles common SCF, Berny, frequency and orientation records, but not every Gaussian link/multilayer ONIOM or Link1 workflow.
6. Scientific validity must be checked using original Gaussian output; synthetic mock energies/frequencies have no physical meaning.
7. The job monitor is single-session and does not recover ownership of a calculation after restarting the GUI. Closing the GUI does not kill owned Gaussian processes.

**Licensing:** GauMon is a separate monitoring utility and does not distribute Gaussian 16 source, binaries, or licensed functionality. The source materials supplied for the original GauMon were described as MIT-licensed; confirm project authorship and redistribution terms before public release. The mock logs are authored for testing only.

## Visual layout acceptance test on Rocky Linux / XRDP

Run this after installing `gaumon-cpu`, to generate screenshots at the intended desktop sizes and perform obvious overlap / clipping checks:

```bash
conda activate gaumon-cpu
python scripts/gui_layout_check.py --output ./qa_screenshots
# For headless CI only:
python scripts/gui_layout_check.py --offscreen --output ./qa_offscreen
```

The command saves **Overview** and **Timeline & ETA** screenshots at 1366×768 and 1920×1080, and a JSON report of detectable card/plot overlap or clipping. The results need to be visually inspected because no automated test can prove that a UI is beautiful or fully accessible across desktop themes and font settings. A Qt graphical runtime was unavailable in the packaging environment, so screenshots were not produced there.

Further details: see [DEEP_SCF_GUIDE.md](DEEP_SCF_GUIDE.md).

## QCNR / LinEq1 vs vibrational frequencies

`LinEq1: Iter=... NonCon=... RMS=... Max=... NDo=...` comes from
Gaussian's iterative linear-equation solve during quadratic SCF, often after
a difficult regular DIIS cycle. It does **not** report vibrations in cm⁻¹.

```bash
python run.py --qcnr-report /path/to/metal_cluster.log
python run.py --freq-report /path/to/optfreq.log
python run.py --qcnr-report examples/stress_logs/qcnr_lineq1.log
python run.py --freq-report examples/stress_logs/harmonic_modes.log
```

The uploaded Fe₂S₂ single-point log exercises L508 and LinEq1 but contains
no `Frequencies --` lines. It should therefore populate Deep QCNR but not
Deep Frequencies. The original uploaded Gaussian logs are not included in the
package; the two new stress fixtures are synthetic.

### Source release validation

See [RELEASE_VALIDATION.md](RELEASE_VALIDATION.md) for test counts, real-log statistics, and remaining workstation-only Qt checks.

### v7.1.1: saved-session startup failure resolved

If the Session Jobs page previously raised `AttributeError: 'NoneType' object has no attribute 'termination'`, update to v7.1.1. The initial draw now treats persisted records as **saved, but not attached**, and only reports completed/failed status when a real log snapshot exists. Opening the GUI never resumes an external Gaussian calculation automatically.

Run the Qt smoke tests on the actual Rocky Linux GUI environment:

```bash
QT_QPA_PLATFORM=offscreen python -m unittest discover -s tests -p test_gui_optional.py -v
./run.sh --gui
```
