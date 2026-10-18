# GauMon Desktop — validation record

The updated source was built in a Python 3.13 container with NumPy and psutil.

- Full automated test suite: **57 evaluated, 55 passed, 2 Qt-only smoke tests skipped** because PyQt5/pyqtgraph were not available.
- Python modules and scripts passed syntax compilation.
- Synthetic fixture `qcnr_lineq1.log`: two quadratic SCF linear-equation solves, the first converged and second unfinished; parser handles chunk boundaries.
- Synthetic fixture `harmonic_modes.log`: six vibrational modes, one imaginary mode, printed reduced masses, force constants, and IR intensities.
- User's original `07_advanced_Fe2S2_toy_SP (copy 1).log.txt`: 8 QCNR orbital iterations, 5 separate LinEq1 solves, 4 converged / 1 incomplete, no harmonic frequencies.
- User's prior `WaterDrops.log` and `d.log` have no QCNR/LinEq1 solver sections in those uploaded snapshots.
- Existing GauMon tests for regular SCF, running/attached log handling, stop semantics, and diagnostics remain passing.
- The proprietary Gaussian output files were used only for local development validation and are **not included in the release ZIP**.

**Not verified:** actual Qt GUI rendering / font clipping, native launcher installation on Rocky Linux, nor Gaussian quantum-chemistry correctness. Screenshot validation script now covers Overview, Timeline, Deep SCF, Deep QCNR, and Deep Frequencies at 1366×768 and 1920×1080. See the quick-install guide.
