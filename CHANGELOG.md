## 7.1.1 — Session Jobs GUI startup hotfix (2026-10-19)

- Fixed immediate startup failure in `gaumon_ext/gui.py` when restored session entries had no attached reader (`AttributeError: 'NoneType' object has no attribute 'termination'`).
- Added `describe_observed_job()` with explicit `SAVED · RESUME TO ATTACH` and `ATTACHING · WAITING FOR LOG` statuses, plus error/unverified/completed states. A saved record is not proof of a running Gaussian process.
- Added five headless regression cases and one Qt startup smoke test for persisted-but-unattached logs. No automatic job or process restart.
- Preserved v7.1 full-screen startup, launch CWD, input-editor autoload, Deep SCF/QCNR/Frequencies, time analytics and safe owned-process control.
- Verified 80 tests from the suite excluding mock-subprocess tests; 76 passed, 4 Qt-only tests skipped (Qt unavailable in the build container). Five focused new status tests pass. Full suite did not complete in this container because an existing mock subprocess test intermittently hangs; it must be repeated on Rocky Linux.

## 7.1.0 — Full-screen / CWD / input intelligence / deque fix

- Startup uses true full-screen display by default. F11 toggles immersive mode; Esc restores the maximized window. Set `GAUMON_START_MAXIMIZED=1` to launch maximized instead, useful on certain XRDP window managers.
- All common file and export dialogs use the current working directory as their default. The shell launcher preserves the directory it was invoked from.
- Starting a new Gaussian job immediately shows its explicit `.gjf`/`.com` file in the Input Checker / Editor. Attaching a log automatically loads an existing same-stem `.gjf`/`.com` beside it; session records preserve explicitly associated original input paths for read-only resume.
- If an original input cannot be verified, the app reports this without silently loading an unrelated file or misrepresenting a route reconstructed from the log as the exact original.
- A dedicated Keyword explanations tab lists recognized Link0 and route keywords, purposes and caveats. Known route tokens appear in green with hover/cursor explanations. Editing updates checks after a short debounce. Existing line numbers, visible whitespace, and Save As remain.
- Fixes the observed recurring `TypeError: sequence index must be integer, not 'slice'` during advanced ETA refresh. Adjacent observations are iterated with `itertools.pairwise`, compatible with `deque`, and regression-tested.
- Automatic input loading never discards unsaved edits without confirmation. No automatic Gaussian job execution, system-wide process kill or daemon discovery is introduced.
- PyQt5 rendering must still be verified on the user's actual Rocky Linux/XRDP environment.

## 7.0.0 — Consolidated GauMon Desktop release

- Rebuilt on the actual uploaded v5.2 ZIP; retained Deep SCF, Deep QCNR and Deep Frequencies.
- Integrated five original v6.1 utility pages into the same QStackedWidget native Qt desktop: Time Graphs, Input Checker/Editor, Advanced ETA, Warnings & Errors, Session Jobs.
- Preserved measured Gaussian `Leave Link ... cpu: ... elap: ...` timing separately from inferences, with CSV scientific-data views unchanged.
- Added persistent, explicit manual resume of session-owned log paths, without automatic external job execution.
- Improved `--env-check` to import `PyQt5.sip` and related GUI modules, pinned `pyqt5-sip` in `environment.yml`.
- Fixed `run.sh` to preserve the caller's working directory; changed new-job and browse defaults to cwd.
- Kept safe stop limited to GauMon-owned Gaussian process groups and protected the main monitor from switching while such a job is active.
- Added headless timing/warnings/input CLI options and consolidated tests for launch-from-anywhere and resumption safety.
- Widened sidebar and compacted its navigation buttons for fifteen-view layout (visual QA on target pending).

## 5.2.0 — Deep QCNR / LinEq1 + Deep Frequencies

- Added optional CPU-only Deep QCNR / LinEq1 workspace for L508 quadratic SCF residuals and outer orbital-iteration history; live inner plots reset per new LinEq1 solve and retain historical solves.
- Added separate optional Deep Frequencies workspace for actual Gaussian `Frequencies --` normal-mode quantities with mode table and CSV export.
- Added stable per-line parsers for D-notation QCNR convergence limits, LinEq1 NonCon/RMS/Max/NDo, frequency values, reduced masses, force constants, and IR/Raman properties.
- Added `--qcnr-report` and `--freq-report` headless commands, plus two synthetic stress log fixtures and regression tests.
- Confirmed the user-supplied Fe₂S₂ SP log has five LinEq1 solves, four converged and one still running, but no harmonic-frequency records.
- GUI window title remains unversioned; no GPU or browser dependencies were added. Extended screenshot QA for both new workspaces.
- 57 tests evaluated: 55 passed, 2 Qt-only smoke tests skipped due to missing PyQt5/pyqtgraph in this build environment.

## 5.1.0 — Deep SCF diagnostics

- Added optional Deep SCF page; hidden until the user requests it, with automatic per-SCF reset and preserved prior-block selector.
- Added streaming parsing of Gaussian `Cycle N Pass N IDiag N` records and DIIS error, Delta-E, RMSU, RMSDP, MaxDP, OVMax, gaps, CP, NSaved, coefficients, and related diagnostic fields.
- Added four CPU-rendered scientific plots, per-iteration table, details panel, and CSV export.
- Fixed duplicated legacy optimizer event ingestion; SCF iteration and optimizer step histories remain separate.
- Added 3 more advanced Gaussian input templates and 4 deterministic synthetic SCF stress logs with a generator and regression tests.
- Added `--deep-report` and `--stress-case` headless CLI modes, optional stress fixture picker in Test Cases.
- Updated native screenshot QA to include deep SCF at 1366×768 and 1920×1080.
- Verified parser against the user's `d.log` and `WaterDrops.log` without bundling original logs.

# Release notes

## Current build

- Refined warm-light, CPU-only native desktop interface with no version badge or logo.
- Added a dedicated **Timeline & ETA** page with four readable forecast cards, an optimization-cycle table, and Gaussian link history.
- Restored event-by-event link and optimizer-step history when attaching to previously running logs, including repeated links.
- Fixed an output-refresh condition that previously depended on the length of an already-full 130-line tail.
- Separated real Gaussian log timestamps from GauMon observation intervals, explicitly declining to estimate elapsed cycles from instantly parsed historical events.
- Added a priority diagnosis for `FormBX had a problem`, preserving the failing atom-number triplets and quartets and linking to Zhe Wang's guide.
- Added a visible failure banner linking to detailed diagnostics without expanding the Overview layout.
- Added `--attach`, `--inspect`, `--list-jobs`, and a user-local installer for the `gaumon` command and desktop application shortcut.
- Kept `gaumon-cpu` as the dedicated Python/Qt environment; the graphical interface remains software-rendered.
