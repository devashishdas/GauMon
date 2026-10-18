# GauMon | Deep SCF diagnostics

The **Deep SCF** workspace is hidden until explicitly opened from Overview or Analysis. It does not affect the Gaussian calculation or the active output file.

## What constitutes one graph block?

An **SCF block** is a sequence beginning at a Gaussian `Cycle N Pass P IDiag I:` line. Each iteration belongs to the same electronic self-consistent-field procedure at a molecular geometry. The block is marked completed only after `SCF Done:`. When Gaussian begins another SCF sequence at `Cycle 1`, the *Follow newest* graphs restart at their first iteration. The previous completed block remains in the dropdown. Gaussian may perform more than one SCF procedure at one geometry-optimization step, and the `Step number` record may appear **after** `SCF Done`; labels are assigned only when supported by the log.

No plotting occurs for the deep page when it is closed. The parser retains at most 120 SCF procedures × 400 iterations each. The page deliberately uses only CPU pyqtgraph/Qt QPainter, not OpenGL or a browser.

## Four charts

| Chart | Data | Display | Interpretation |
|---|---|---|---|
| SCF electronic energy | `E=` within each cycle | Hartree, linear axis | Electronic energy trend, not optimization RMS force |
| Absolute energy difference | `Delta-E=`; falls back to `DE=` | Hartree, positive magnitude on log axis | Energy stabilization; signs preserved in table |
| DIIS / RMSU | `DIIS: error=`, `RMSU=` | Two independent diagnostic series, log axis | Iterative convergence; neither is directly a bond displacement |
| Density update | `RMSDP=`, `MaxDP=` | Positive magnitude on log axis | Changes in the density matrix |

Zero or missing values are omitted from logarithmic charts. Plot colors are descriptive, not a Gaussian convergence pass/fail verdict.

## More fields available by selecting a table row

`Pass`, `IDiag`, `NSaved`, `CP`, `Coeff-Com`, `Coeff-En`, `Coeff`, `Gap`, `WtCom`, `WtEn`, `ErrMin`, `ErrMax`, `BMatC`, `BMatP`, `EMaxC`, `OVMax`, `Rises`, and `Damp` are retained when present. Gaussian can place adjacent fixed-width numbers on one line (e.g. `0.274D-01-0.804D-01`); the parser uses number-token recognition rather than whitespace splitting to handle this.

## Fast verification

```bash
# Uses Python only; no GUI required
python run.py --deep-report /path/to/d.log
python run.py --deep-report /path/to/WaterDrops.log

# Synthetic examples that ship with the package
python run.py --stress-case many_steps
python run.py --stress-case oscillation
python run.py --stress-case incomplete
python run.py --stress-case chunked

# Opens Qt GUI and simulates Gaussian-like cycle output
python run.py --demo-ui
# Then click: Overview > Deep SCF diagnostics
```

Results observed from the development samples: `d.log` contains one unfinished 19-cycle SCF procedure, and `WaterDrops.log` contains 19 completed SCF procedures, 18 with explicitly printed optimizer-step labels. These counts are properties of the uploaded snapshots, not proof that the calculations converged or finished. They do not measure elapsed SCF time because Gaussian does not timestamp every internal iteration.

## Harder scientific test inputs

The **Test cases** panel lists four simple educational inputs plus `trimer`, `stretched`, and `fe2s2` advanced inputs. Only `--run-real` or a deliberate GUI launch starts licensed Gaussian. The Fe–S example is an illustrative spin-sensitive electronic-structure stress case, **not** a validated broken-symmetry calculation. The seven templates are input-syntax examples, not guaranteed successful chemistry or optimized structures.

## Limits and safety

- An unfinished Gaussian log may end in the middle of a cycle. Such a cycle can have missing metrics; the UI must not invent values or mark the SCF as complete.
- `SCF Done` is not equivalent to *geometry optimization completed*; the two statuses are separate.
- The deep plot does not itself calculate SCF forces, total CPU usage, a reaction barrier, or real physical convergence criteria.
- This build does not contain the uploaded proprietary Gaussian output files. Its synthetic fixtures were generated separately for testing.
- Current graphical layout and font clipping must be verified on the actual PyQt5 workstation using `python scripts/gui_layout_check.py --output ./qa_screenshots`.
