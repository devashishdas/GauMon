# GauMon Desktop v7.1.1 — install and run

## 1. Unpack in a permanent location

```bash
unzip GauMon_Desktop_CPU_v7_1_1.zip
cd GauMon_Desktop_CPU_v7_1_1
```

Avoid extracting over an older v5.2 or v6 folder. Existing files and Gaussian
jobs are not modified during setup.

## 2. Update existing CPU environment (includes PyQt5 SIP)

```bash
conda env update -n gaumon-cpu -f environment.yml --prune
conda activate gaumon-cpu
python run.py --env-check
```

If the environment does not yet exist, use
`conda env create -f environment.yml` instead.

In case the checker reports `No module named 'PyQt5.sip'`, run:

```bash
conda install -n gaumon-cpu -c conda-forge \
  --override-channels --strict-channel-priority \
  --force-reinstall pyqt pyqt5-sip pyqtgraph
python -s run.py --env-check
```

If imports still fail, check for pip/Conda conflicts; do not manually copy
`.so` files between Python versions.

## 3. Test and launch

```bash
python run.py --test
./run.sh
```

The default working directory stays where you invoked `./run.sh`; it no longer
changes to the package folder. Browse for any existing `.gjf` or `.log`.

## 4. Install command to run anywhere

```bash
./install_user.sh
export PATH="$HOME/.local/bin:$PATH"
gaumon --gui
# From inside a directory with calculation.log:
gaumon --attach "$PWD/calculation.log"
```

This replaces a recognized previous GauMon per-user launcher only after
checking the interpreter and both ownership markers. It doesn't use sudo.

## 5. Useful headless commands

```bash
python run.py --timing-report /path/to/job.log
python run.py --warnings-report /path/to/job.log
python run.py --check-input /path/to/job.gjf
python run.py --deep-report /path/to/job.log
python run.py --qcnr-report /path/to/job.log
python run.py --freq-report /path/to/job.log
python scripts/gui_layout_check.py --output ./qa_screenshots
```

The GUI's **Session jobs** page stores only the paths you explicitly opened or
launched. After restarting GauMon choose **Resume selected** to reopen those
logs read-only; it never restarts Gaussian itself.

## Start fullscreen, use current directory and automatically inspect inputs

Start GauMon **from your calculation directory**. That directory is used for
new jobs, opening inputs/logs and CSV exports by default:

```bash
cd /home/exampleuser/my_calculations
gaumon --gui
# Or directly from the extraction folder, while preserving the current dir:
/path/to/GauMon_Desktop_CPU_v7_1_1/run.sh --gui
```

The window starts in **true full-screen** mode. Press **F11** to toggle
full-screen or **Esc** to restore a maximized window. For certain XRDP window
managers use `GAUMON_START_MAXIMIZED=1 ./run.sh --gui` instead.

When you **start a new job**, GauMon loads the exact input you selected in
**Input checker / editor**. When you **attach a log**, it tries an existing
same-named `.gjf`/`.com` next to that log. When no matching original exists,
it explains why the source cannot be reliably inferred from the output; use
**Open .gjf**. Session records preserve any explicit input association for
later read-only resumption.

The Input Checker includes syntax highlighting, a live keyword reference tab,
visible whitespace, line numbers, structural warnings and an explicit Save As.
Recognized keywords such as `#p`, `opt`, `freq`, `scf=xqc`, `maxcycle`,
`%nprocshared` and `%mem` are explained in plain English. The checker is
not a substitute for Gaussian's chemistry and input validation.

## Regression verification

```bash
python run.py --test
python -m unittest -v tests.test_input_context_v71
python scripts/gui_layout_check.py --output ./qa_screenshots
```
