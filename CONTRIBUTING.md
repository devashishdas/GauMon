# Contributing to GauMon

GauMon is an MIT-licensed desktop tool for Gaussian output monitoring. Issues, bug reports, documentation updates, and reproducible tests are welcome.

## Development workflow

1. Fork the repository and create a feature branch. Keep changes limited and explain the behavior they modify.
2. Install either `environment.yml` (Conda) or `requirements.txt` into a private venv. Do not run package installation using `sudo pip`.
3. Run `python run.py --test` and `python -m unittest discover -s tests -p 'test_install_wizard.py' -v`.
4. If editing the GUI, test `QT_QPA_PLATFORM=offscreen python scripts/gui_layout_check.py --offscreen --output ./qa` and then inspect screenshots using a real desktop at 1366×768 and 1920×1080.
5. Submit a pull request describing the scientific parsing behavior, tests, any changed dependencies, and relevant screenshots.

Never attach proprietary Gaussian binaries, license files, private research logs, or unredacted cluster credentials. Use synthetic or explicitly redistributable fixtures. Gaussian jobs must never be launched during automatic tests.

## Bug reports

Include GauMon version, Linux distribution and version, Python/Conda version, install mode, the exact terminal command, `gaumon --env-check` output, and the relevant exception. Redact input files, usernames, paths, molecular structures, and Gaussian output unless you have permission to publish them.

## Support

Use [GitHub issues](https://github.com/devashishdas/GauMon/issues) for usage questions and reproducible defects. Pull requests are reviewed for test coverage, safety, maintainability, and source-license compatibility.
