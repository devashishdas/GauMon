# Deep QCNR / LinEq1 and Deep Frequencies

## The scientific distinction

`LinEq1` is **not a frequency calculation**. In Gaussian 16 it appears during the
quadratically convergent SCF process (e.g. L508 when `SCF=XQC` switches over after
ordinary SCF struggles). It reports inner iterations of a **linear-equation
solver** for the orbital optimization step. These nested solver iterations are
neither geometry-optimization steps nor vibrational normal modes.

Typical output:

```
(Enter .../l508.exe)
Quadratic Convergence SCF Method.
Iteration   1 EE= -4918.03385767531 Grad=9.105D-03 Max rot= 0.0996 deg.
QCNR: CnvC1=1.32D-08 CnvC2=1.32D-07
LinEq1: Iter= 0 NonCon= 1 RMS=2.71D-04 Max=1.32D-02 NDo=1
...
LinEq1: Iter=35 NonCon= 0 RMS=1.82D-09 Max=4.63D-08 NDo=1
Linear equations converged to 1.319D-08 1.319D-07 after 35 iterations.
```

Meanings: `Iter` is the inner solver iteration counter; `NonCon` is the number
of unconverged RHS/equations; `RMS` and `Max` are linear-equation residual
magnitudes; `NDo` is the number of linear systems to process in that printout.
`QCNR` outer `Iteration` contains electronic energy (`EE`), energy change,
orbital gradient (`Grad`), and maximum rotation in **degrees**. The gradient
and rotation are not plotted on a shared numerical axis (different units).
A solve with `NonCon=0` is *converged at the inner-linear-solver level*; it does
**not** establish total SCF convergence. Log EOF before convergence is *running
or incomplete* rather than successful.

## Deep QCNR view

On Overview click **Deep QCNR / LinEq1**. This optional view has four charts:

- LinEq1 RMS residual versus its **own** iteration (logarithmic Y, positive only)
- LinEq1 maximum residual versus its own iteration (logarithmic Y)
- Orbital-optimization EE versus QCNR outer iteration (Hartree)
- Orbital-gradient magnitude versus QCNR outer iteration (printed units)

Select prior linear solve runs; the live residual plots reset when a new
LinEq1 solve starts. Table includes `Iter`, `NonCon`, `RMS`, `Max`, `NDo`.
The readout includes requested convergence limits and max orbital rotation.
You can export the selected solve to CSV.

## Deep Frequencies view

On Overview click **Deep Frequencies**. This genuinely separate analysis only
uses Gaussian `Frequencies --` lines and their associated `Red. masses --`,
`Frc consts --`, `IR Inten --`, and optional Raman/depoloarization records.
It shows (i) normal-mode frequency against mode index and (ii) IR intensity
against frequency where intensities are present. The IR chart is **not** a
broadened/normalized experimental spectrum: intensities are unmodified printed
numbers. Modes with negative frequencies are counted as imaginary.

An SP/SCF calculation such as the supplied Fe2S2 example has **no vibrational
frequency modes**, even when LinEq1 output is present.

## Test commands (no licensed Gaussian invocation)

```bash
conda activate gaumon-cpu
python run.py --test
python run.py --qcnr-report examples/stress_logs/qcnr_lineq1.log
python run.py --freq-report examples/stress_logs/harmonic_modes.log
python run.py --qcnr-report '/path/to/07_advanced_Fe2S2_toy_SP (copy 1).log.txt'
python run.py --demo-ui
```

Open a log directly from anywhere (after running `./install_user.sh`):

```bash
gaumon --attach /path/to/calculation.log
```

## Validation with the supplied Fe–S example

The particular uploaded 3263-line snapshot contains eight QCNR outer orbital
iterations; five separate LinEq1 solves: **four completed inner solves and a
fifth truncated during its 16th printed iteration** (iteration indices 0–15).
There are **zero harmonic `Frequencies --` lines**. The plain DIIS stage says
its convergence criterion was not met before Gaussian entered quadratic SCF.
These facts do not by themselves mean the full job subsequently finished.

No proprietary Gaussian log, executable, or checkpoint is distributed in
this source package. The included stress logs are independently authored
synthetic fixtures. Neither graphical page needs a GPU.

## GUI verification

Run `python scripts/gui_layout_check.py --output ./qa_screenshots` on your
Rocky Linux machine with PyQt5 in `gaumon-cpu`. Check the Deep QCNR and Deep
Frequencies screenshots at 1366×768 and 1920×1080. Headless code tests
cannot verify actual Qt text clipping or visual layout on your workstation.
