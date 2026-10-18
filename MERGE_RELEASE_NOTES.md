# GauMon Desktop v7 — merged feature release

GauMon Desktop v7 is a **single source distribution based on the user's original
v5.2 ZIP**, with the requested v6.1 utility pages integrated into the existing
native v5.2 `MainWindow`. It deliberately avoids embedding the v6.1 legacy
`GauMon.py` front end or a second whole desktop application.

## Fifteen pages in one interface

### Existing v5.2 pages, preserved

1. Overview: six live metrics and four scientific plots.
2. Run calculation: authorized `g16` execution, start/stop (managed processes
   only), read-only log attachment and PID selection.
3. Timeline & ETA: SCF optimization steps and link history.
4. Analysis: SCF table, frequency and geometry diagnostics.
5. Geometry: 2D projection with Kabsch-based geometry metrics.
6. Live log: log tail and source-credited FormBX/error diagnosis.
7. Test cases: seven simple/advanced Gaussian inputs and synthetic jobs.
8. Deep SCF: on-demand electronic-cycle plots, DIIS/RMS residuals, CSV.
9. Deep QCNR: on-demand L508/LinEq1 diagnostics, CSV.
10. Deep Frequencies: actual normal-mode frequencies/intensities, CSV.

### New integrated utility pages

11. **Time graphs:** measured Gaussian **Leave Link** elapsed and CPU seconds,
    link aggregates, sequence and per-link contributions. Work still executing
    within an unfinished link is *not* assigned fabricated runtime; these are
    computation-stage, not per-atom, timings. Plain-English interpretation is
    displayed above and below the graphs.
12. **Input checker & editor:** .gjf/.com opener, colored route/Link0/Link1 and
    charge/multiplicity records, monospaced editing, numbered gutter,
    configurable visible whitespace and line endings, error-line highlights,
    clickable problems table, Save As. Heuristic check, not a Gaussian compiler.
13. **Advanced time estimator:** actual newly observed optimization-step timing,
    near-term +10/50/100/150-step forecasts and descriptive pace ranges.
    **Backfilled logs alone do not justify historical step timing.** No
    calibrated uncertainty interval or certain job-finish promise is claimed.
14. **Warnings & errors:** structured log-line evidence, FormBX, SCF failure,
    normal/error termination and conservative suggested checks.
15. **Session jobs:** explicitly enrolled job logs, up to 20 active readers,
    local JSON session metadata and manual **Resume** after reopening GauMon.
    No automatic external calculation start and no process-wide kill.

## Important behavior and limits

- The v5.2 monitor remains the primary live dashboard. The new analytics
  pages use a separate *read-only incremental timing parser* on explicitly
  enrolled paths. They can take several refreshes to catch up on large logs.
- Switching the primary monitor away from a running GauMon-owned job is
  blocked, retaining safe process ownership controls. Reading a saved log
  never assumes the external Gaussian process is alive.
- A saved job path is not a saved Gaussian process or checkpoint. Manual Resume
  reopens its existing `.log`; it does not resume the calculation numerically.
- The global process `pkill` behavior from the original GauMon.py is not
  part of this release; stop actions only target a GauMon-owned process group.
- Additional `--timing-report`, `--warnings-report`, and `--check-input`
  commands run without Qt or a Gaussian executable.
- Native Qt may need a clean Conda SIP installation. `--env-check` actually
  imports `PyQt5.sip` and reports its path instead of relying on `find_spec`.
- `run.sh` preserves the **calling directory**. File-dialog fallbacks and new
  job working directory use the user's current working directory.

## Real-log sanity checks

The attached developer examples were read locally, **not bundled** in v7:

- `WaterDrops.log`: 19 completed L502 link calls, 93.80 s accumulated
  elapsed time in L502; entire log contains 212 completed link records.
- `07_advanced_Fe2S2_toy_SP (copy 1).log.txt`: L502 spent 152.30 s elapsed;
  the interrupted/in-progress L508 stage has no reported completed-link
  duration in the supplied snapshot.

**No physical atom-by-atom wall-time attribution** is claimed. Electronic
structure calculations couple atomic orbitals across the system.

## What's not yet independently validated

The build environment has no PyQt5/pyqtgraph, so actual Qt screenshots and
click-throughs at 1366x768 / 1920x1080 remain unverified. Run:

```bash
python scripts/gui_layout_check.py --output ./qa_screenshots
```

This QA script saves screenshots of Overview, Timeline, the three deep-analysis
pages, and the five new utilities. Review them manually for text clipping and
window behavior under Rocky Linux / XRDP. Also verify launcher installation,
editor save and real-world large-log performance on the target workstation.
