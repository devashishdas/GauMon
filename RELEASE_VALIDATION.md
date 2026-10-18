# GauMon Desktop v7 — source and test verification

## Verified in this build

- **78 unit/integration tests executed; 76 passed, 2 skipped.** The two
  skips are Qt desktop smoke tests because PyQt5 and pyqtgraph are not installed
  in this Python builder environment. Test run is recorded as a local console
  test log; the package can rerun it using `python run.py --test`.
- Python modules and scripts compile, and launch/installer shell scripts pass
  `bash -n` syntax verification.
- CLI starts from another current working directory using `run.sh` and a
  **relative** log name. It does not change its caller's working directory.
- `--env-check` attempts to import `PyQt5.sip`, `PyQt5.QtCore`, `PyQt5.QtGui`,
  `PyQt5.QtWidgets`, `pyqtgraph`, NumPy and psutil; missing SIP is now diagnosed
  explicitly with conda-forge recovery commands.
- Gaussian stage elapsed/CPU totals are measured from Gaussian's reported
  complete `Leave Link` records and not attributed to individual atoms.
- Input checker detects route/title separator errors and shows specific lines,
  without running Gaussian.
- Session store persists only manually enrolled log paths and never contains
  executable recovery commands or automatic job restarts.
- Existing Deep SCF, QCNR/LinEq1, harmonic frequency, FormBX/L103 diagnosis,
  and synthetic scientific tests remain included.
- New tests cover integrated CLI reporting, session persistence, cwd safety,
  Qt/SIP dependency diagnostics and absence of system-wide Gaussian killing.

## User log sanity checks (not redistributed)

- `WaterDrops.log`: L502 19 calls, 93.80 seconds elapsed in completed L502
  operations; 212 complete link records across the supplied snapshot.
- `07_advanced_Fe2S2_toy_SP (copy 1).log.txt`: L502 152.30 seconds elapsed
  in one completed call; no fabricated completed L508 time.

## What is *not* verified

- **Qt GUI appearance, click-through behavior, and integration with the
  window manager** remain untested on Rocky Linux 8 / XRDP. This is a source
  release pending that test, not a verified production desktop build.
- End-to-end integration with a real licensed Gaussian process is not tested
  in this container. It is not bundled or executed by the test suite.
- ETA ranges are uncalibrated; the app makes no guaranteed completion-time
  prediction or per-atom CPU-time attribution.
- Very large logs and 20 simultaneously attached logs need workstation-level
  performance profiling, particularly during initial backfill.

For visual checks after installation:

```bash
python run.py --env-check
python scripts/gui_layout_check.py --output ./qa_screenshots
```

The screenshot generator covers Overview, Timeline, Deep SCF, Deep QCNR,
Deep Frequencies, Time Graphs, Input Checker, Advanced ETA, Warnings & Errors,
and Session Jobs at two display resolutions. Review any clipped labels and
check multiple log sessions manually before publishing a public release.
