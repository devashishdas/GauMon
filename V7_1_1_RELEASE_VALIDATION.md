# GauMon v7.1.1 GUI startup hotfix

## Source of failure

User's real Rocky Linux log showed `JobsPage.redraw()` on session reload at
`gaumon_ext/gui.py` line 439 attempting `st.termination` where `st` was `None`.
This happens when a saved session is loaded but its log has **not** been reattached.
Because the exception occurs in `MainWindow.__init__()`, it prevents the GUI
opening and causes unrelated Qt smoke tests to fail. It does not indicate a
Gaussian error or missing Qt libraries.

## Fix

`gaumon_ext/session.py` now exports `describe_observed_job(engine)`.
`gaumon_ext/gui.py` uses this for each row and guards all snapshot access.
Saved sessions show `SAVED · RESUME TO ATTACH`. A newly-created reader
without a snapshot displays `ATTACHING · WAITING FOR LOG`. Other statuses
are based on observed errors, terminations, or log growth, with process
verification explicitly withheld. This maintains read-only manual resume.

## Validation performed in container

- All Python source and tests compiled successfully.
- Five focused session-status regression tests passed.
- 80 tests from the whole suite excluding mock-subprocess TestMockCases:
  76 passed, 4 Qt-only skipped because PyQt5/pyqtgraph are not installed.
- Full suite was attempted twice but intermittently stalled in unrelated
  existing mock subprocess cases; no comprehensive pass is claimed.
- Python GUI source could not be executed visually in this environment.
- Clean extraction and focused tests will be repeated on the final ZIP.

## Verify on the user's workstation

```bash
conda activate gaumon-cpu
cd GauMon_Desktop_CPU_v7_1_1
python run.py --env-check
QT_QPA_PLATFORM=offscreen python -m unittest discover -s tests -p test_gui_optional.py -v
./run.sh --gui
```

Using `QT_QPA_PLATFORM=offscreen` is only for the Qt tests, not for the actual
interactive window. The project still opens full-screen by default, with F11
or Esc to change view modes. This package never launches real Gaussian
calculations without explicit opt-in.
