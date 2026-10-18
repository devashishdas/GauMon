"""Headless Qt dependency preflight for the GauMon desktop application."""
from __future__ import annotations
import importlib
import os
import site
import sys
import traceback

REQUIRED = ('PyQt5.sip', 'PyQt5.QtCore', 'PyQt5.QtGui', 'PyQt5.QtWidgets', 'pyqtgraph', 'numpy', 'psutil')


def qt_environment_report(*, verbose=True):
    """Actually import GUI packages; a find_spec('PyQt5') test is insufficient."""
    lines = [f'Python executable: {sys.executable}',
             f'Python version: {sys.version.split()[0]}',
             f'CONDA_PREFIX: {os.environ.get("CONDA_PREFIX") or "(not set)"}',
             f'PYTHONPATH: {os.environ.get("PYTHONPATH") or "(not set)"}',
             f'User site enabled: {site.ENABLE_USER_SITE}',
             f'User site directory: {site.getusersitepackages()}']
    ok = True
    for name in REQUIRED:
        try:
            module = importlib.import_module(name)
            location = getattr(module, '__file__', None) or '(namespace / built-in)'
            lines.append(f'OK      {name}: {location}')
            if name == 'PyQt5.QtCore':
                lines.append('        Qt: ' + module.QT_VERSION_STR + '  PyQt: ' + module.PYQT_VERSION_STR)
        except Exception as exc:
            ok = False
            lines.append(f'FAILED  {name}: {type(exc).__name__}: {exc}')
            if verbose:
                last=''.join(traceback.format_exception_only(type(exc),exc)).strip()
                if last not in lines[-1]:lines.append('        '+last)
    if not ok:
        lines.extend([
            '', 'Suggested conda-forge repair (inside gaumon-cpu):',
            '  conda install -n gaumon-cpu -c conda-forge --override-channels --strict-channel-priority --force-reinstall pyqt pyqt5-sip pyqtgraph',
            'Check mixed pip/conda packages with:',
            '  python -m pip show PyQt5 PyQt5-sip PyQt5-Qt5',
            'Test for user-site interference with:',
            '  python -s run.py --env-check',
            'Avoid manually copying or symlinking sip extensions across different Python versions.',
        ])
    return ok,'\n'.join(lines)
