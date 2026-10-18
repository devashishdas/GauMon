# GauMon Desktop CPU v7.1 — release validation

## Fixed with evidence

- Corrected the Advanced ETA `TypeError: sequence index must be integer, not 'slice'`: an underlying `collections.deque` was sliced via `pairs[1:]` in a Qt redraw callback. Timing sequences are now handled with `itertools.pairwise` in a headless, tested helper (`gaumon_ext.analytics.adjacent_step_durations`). The original GUI exception could repeat during periodic redraw without failing the non-GUI test suite.
- Default GUI display now uses `showFullScreen()`. F11 toggles; Esc restores maximized mode. Users can set `GAUMON_START_MAXIMIZED=1` for a maximized initial window instead (XRDP fallback).
- Current launch working directory is the default new-job directory and the default base for input, log, exported examples, and CSV file choosers. `run.sh` and installed per-user launcher preserve the caller's working directory.
- Starting a Gaussian job through GauMon automatically loads the user-selected `.gjf` or `.com` into the input editor. Attaching or selecting a session log attempts an existing same-stem `.gjf` or `.com` beside the output. Saved session records retain an explicitly associated original path for later manual read-only resume.
- Missing input is explicitly disclosed rather than fabricating input from a Gaussian log. Unsaved input edits must be approved before automatic replacement; the editor uses Save As.
- Colored route keywords, `%` link-zero keywords, table of plain-English keyword explanations, a contextual cursor-based guide, and debounced validation augment existing colored input syntax, line numbers, visible spaces, tabs and line endings.

## Automated checks

- `python run.py --test` validates GauMon parser, mock-process safety, existing deep scientific pages, event-time analytics and the updated editor/association code where Qt is available.
- `python -m unittest -v tests.test_input_context_v71` independently tests the deque regression, scientific keyword extraction, adjacent `.gjf`/`.com` association, multiple Link1 sections and session-provenance persistence. No Gaussian executable or GPU is needed.
- Python bytecode compilation of core, utility, CLI, scripts and tests, and `bash -n` shell checks pass.
- The attached `WaterDrops.log` and `d.log` each have a matching input in this development environment, and both produce recognized keywords. The uploaded Fe2S2 `*.log.txt` lacks its original `.gjf` beside it; the application must not guess.
- After startup on Rocky Linux run `python scripts/gui_layout_check.py --output ./qa_screenshots` and test the new editor interaction and full-screen keyboard shortcuts. PyQt5/pyqtgraph were unavailable in this build environment, so actual GUI rendering is **not** claimed as locally verified.

## Important scope

There is no automatic Gaussian run on launch or attach; any real Gaussian execution still requires explicit Start. Merely opening an input in the editor is read-only. Session resumes attach to logs read-only and cannot kill unrelated processes. Keyword explanations are an aid, not scientific validation of methods, spin state, basis selection or syntax.
