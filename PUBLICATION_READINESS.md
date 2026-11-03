# GauMon — publication and release readiness

This is a release QA checklist, **not** an assertion that a journal has accepted the software.

## Scientific scope

- Standalone user-space monitor and parser for Gaussian 16 inputs and outputs.
- Live and historical monitoring, SCF/optimization/frequency diagnostics, timing/ETA estimates, read-only attachment, and user-approved process control.
- No licensed Gaussian binary distributed; synthetic logs and mock computations are only test fixtures and must never be presented as physical results.

## Evidence completed in the repository

- MIT LICENSE; README; synthetic and advanced example fixtures.
- Python unit tests for parsing, UI diagnostics, command-line behavior, and sessions.
- A Linux wizard with Conda and private `venv` modes, explicit PATH opt-in, and a stable launcher.
- Conda-forge and pip dependency specifications; GitHub Actions test workflow.
- Citation metadata and contribution/issue-reporting guidance.

## Release gates requiring manual, external verification

- [ ] Test **fresh** Ubuntu and Rocky Linux installations, both Conda and venv, with normal user permissions.
- [ ] Verify no `sudo` or globally installed Python packages are required by the wizard.
- [ ] Test normal desktop and XRDP/X11, `gaumon --gui`, attach, mock-case, layout at 1366×768 and 1920×1080.
- [ ] Audit the tests and record successful CI run IDs and screenshots.
- [ ] Confirm parsing against consented, real Gaussian logs (SP, optimization, frequency, TS, failed, resumed) and document known gaps.
- [ ] Check permissions for all icons, screenshots, and example outputs before including them in publication.
- [ ] Verify multi-user usage, update/reinstallation behavior, release artifacts, and support instructions.
- [ ] Create a tagged GitHub release, archive it in an archival repository such as Zenodo, and update the citation version and DOI.
- [ ] Review documentation accuracy on a system not used for development.

## JOSS-specific work

Check the current [JOSS submission criteria](https://joss.readthedocs.io/en/latest/submitting.html) before submission. JOSS currently requires at least six months of **public** development history and evidence of open development, rather than a newly public one-time upload. The GauMon public repository was created on **2026-10-18**; the six-month anniversary would be **2027-04-18**, but time alone does not establish eligibility. A reviewer also needs demonstrated usage, a maintainable and appropriately packaged codebase, examples, community contribution routes, a paper covering the software's need and state of the field, design tradeoffs, impact, and transparent AI-use disclosure.

**Do not assert originality or a 'first' claim without comparing related Gaussian monitoring and log analysis tools.** Avoid claiming production readiness until the fresh-machine tests above have been recorded. The JOSS paper should describe software and its use, not present new scientific results as the core contribution.

### Linux executable-mode and staging regression (2026-11-03)

- [ ] Commit `examples/mock_g16` with executable mode `100755` using
  `git add --chmod=+x examples/mock_g16`.
- [ ] Confirm `bash run.sh --test` works from a newly cloned checkout.
- [ ] Confirm `bash install.sh --mode conda --yes --bashrc --desktop` completes
  all tests on Rocky Linux without missing `run.sh` or `PermissionError`.
- [ ] Confirm `bash install.sh --mode venv --yes --no-bashrc --no-desktop`
  completes on a fresh machine without any Conda environment.
- [ ] Confirm `~/.local/bin/gaumon --gui` launches with a working display.
