#!/usr/bin/env python3
"""GauMon user installation wizard (Linux, no root privileges required).

Installs a private conda-prefix or venv plus a stable source copy; does not
manage Gaussian or touch Gaussian data, scratch, or configuration.
"""
from __future__ import annotations

import argparse
import os
import stat
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import uuid

APP_MARKER = ".gaumon-managed-app"
COMMAND_MARKER = "# GAUMON_MANAGED_LAUNCHER_V2"
COMMAND_ALIASES = ("Gaumon", "GauMon")
DESKTOP_MARKER = "X-GauMon-Installed=true"
BASH_BEGIN = "# >>> GauMon PATH >>>"
BASH_END = "# <<< GauMon PATH <<<"
SOURCE = Path(__file__).resolve().parent


def log(message: str) -> None:
    print(f"[GauMon installer] {message}", flush=True)


def fail(message: str) -> None:
    raise RuntimeError(message)


def select_conda(explicit: str | None) -> str | None:
    if explicit:
        candidate = shutil.which(explicit) or explicit
        return str(Path(candidate).expanduser()) if Path(candidate).is_file() else None
    candidate = shutil.which("conda")
    if candidate:
        return candidate
    for directory in ("miniconda3", "miniforge3", "mambaforge", "anaconda3", "anaconda", "opt/miniconda3"):
        path = Path.home() / directory / "bin" / "conda"
        if path.is_file():
            return str(path)
    return None


def subprocess_checked(argv: list[str], *, env: dict[str, str] | None = None) -> None:
    log("RUN: " + shlex.join(argv))
    subprocess.run(argv, check=True, env=env)


def ask_bool(prompt: str, default: bool = False) -> bool:
    suffix = "Y/n" if default else "y/N"
    choice = input(f"{prompt} [{suffix}]: ").strip().lower()
    return default if not choice else choice in ("y", "yes")


def pick_mode(mode: str, conda: str | None, interactive: bool) -> str:
    if mode != "auto":
        if mode == "conda" and not conda:
            fail("Conda not found. Install Miniforge, specify --conda-exe, or use --mode venv.")
        return mode
    if not interactive:
        return "conda" if conda else "venv"
    print("\nEnvironment choices:")
    print("  1. Conda (recommended when Conda/Miniforge is available)")
    print("  2. Python venv (without Conda; needs Python >=3.10 and pip)")
    response = input(f"Choose [1/2] (default {'1' if conda else '2'}): ").strip()
    result = "venv" if response == "2" else "conda" if response == "1" else "conda" if conda else "venv"
    if result == "conda" and not conda:
        fail("Conda not found; choose Python venv or install Conda first.")
    return result


def clean_python_env() -> dict[str, str]:
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONHOME", None)
    # A user's active Conda session or third-party Qt app can otherwise inject
    # incompatible Qt plugins into GauMon's dedicated environment.
    for name in ("CONDA_PREFIX", "CONDA_DEFAULT_ENV", "QT_PLUGIN_PATH", "QT_QPA_PLATFORM_PLUGIN_PATH"):
        env.pop(name, None)
    env["PYTHONNOUSERSITE"] = "1"
    return env


def check_python(python: str) -> None:
    subprocess_checked([python, "-c", "import sys; assert sys.version_info >= (3, 10), 'GauMon requires Python 3.10+'"])


def environment_python(mode: str, prefix: Path) -> Path:
    return prefix / ("conda-env" if mode == "conda" else "venv") / "bin" / "python"


def create_environment(mode: str, prefix: Path, stage: Path, conda: str | None, python: str) -> Path:
    interpreter = environment_python(mode, prefix)
    if mode == "conda":
        assert conda is not None
        folder = interpreter.parent.parent
        if interpreter.is_file():
            subprocess_checked([conda, "env", "update", "--prefix", str(folder), "--file", str(stage / "environment.yml"), "--prune"])
        else:
            if folder.exists() and any(folder.iterdir()):
                fail(f"Nonempty unrecognized Conda environment directory: {folder}")
            subprocess_checked([conda, "env", "create", "--prefix", str(folder), "--file", str(stage / "environment.yml")])
    else:
        check_python(python)
        folder = interpreter.parent.parent
        if not interpreter.is_file():
            if folder.exists() and any(folder.iterdir()):
                fail(f"Nonempty unrecognized virtual environment directory: {folder}")
            subprocess_checked([python, "-m", "venv", str(folder)])
        if not interpreter.is_file():
            fail(f"Python virtual environment incomplete: {interpreter}")
        subprocess_checked([str(interpreter), "-m", "pip", "install", "--disable-pip-version-check", "-r", str(stage / "requirements.txt")])
        subprocess_checked([str(interpreter), "-m", "pip", "check"])
    if not interpreter.is_file():
        fail(f"Environment has no Python executable: {interpreter}")
    return interpreter


def check_managed(path: Path, marker: str) -> None:
    if not path.exists() and not path.is_symlink():
        return
    if not path.is_file() or path.is_symlink():
        fail(f"Refusing to overwrite unrelated path: {path}")
    content = path.read_text(encoding="utf-8", errors="replace")
    valid_markers = (marker,)
    # Existing GauMon 7.1.1 installations used this explicit owner marker.
    # Treat it as managed for a safe upgrade; never overwrite generic scripts.
    if marker == COMMAND_MARKER:
        valid_markers += ("# GAUMON_USER_LAUNCHER",)
    if not any(tag in content for tag in valid_markers):
        fail(f"Refusing to overwrite unrelated file: {path}")


def atomic_text(path: Path, content: str, executable: bool = False, mode: int | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=".gaumon-", text=True)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(content)
        os.chmod(tmp, mode if mode is not None else (0o755 if executable else 0o644))
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def shell_command(python: Path, app: Path) -> str:
    # shlex.quote prevents whitespace and shell metacharacter expansion.
    return f"#!/usr/bin/env bash\n{COMMAND_MARKER}\nset -euo pipefail\nunset PYTHONPATH PYTHONHOME CONDA_PREFIX CONDA_DEFAULT_ENV QT_PLUGIN_PATH QT_QPA_PLATFORM_PLUGIN_PATH\nexport PYTHONNOUSERSITE=1\nexec {shlex.quote(str(python))} -s {shlex.quote(str(app / 'run.py'))} \"$@\"\n"


def desktop_entry(launcher: Path) -> str:
    # XDG Desktop Entry Exec quoting of one absolute executable path.
    quoted = str(launcher).replace("\\", "\\\\").replace('"', '\\"').replace("`", "\\`").replace("$", "\\$")
    return ("[Desktop Entry]\nVersion=1.0\nType=Application\nName=GauMon\n"
            "Comment=Gaussian 16 calculation monitor\n"
            f'Exec="{quoted}" --gui\nTryExec={quoted}\n'
            "Terminal=false\nCategories=Science;Education;\nStartupNotify=true\n"
            f"{DESKTOP_MARKER}\n")


def alias_command(launcher: Path) -> str:
    # The short launcher itself is the only place aware of the Python environment.
    return f"#!/usr/bin/env bash\n{COMMAND_MARKER}\nexec {shlex.quote(str(launcher))} \"$@\"\n"


def bashrc_command() -> str:
    return 'export PATH="$HOME/.local/bin:$PATH"  # for this terminal only; persistent setup is opt-in'


def validate_bashrc(bashrc: Path) -> None:
    if bashrc.is_symlink():
        fail(f"Refusing to modify symlinked ~/.bashrc: {bashrc}")
    if bashrc.exists() and not bashrc.is_file():
        fail(f"Refusing to modify non-file ~/.bashrc: {bashrc}")
    existing = bashrc.read_text(encoding="utf-8") if bashrc.exists() else ""
    if (BASH_BEGIN in existing) != (BASH_END in existing):
        fail(f"Incomplete GauMon managed section in {bashrc}; inspect manually")


def update_bashrc(bashrc: Path) -> bool:
    """Append one guarded PATH block, preserving existing bytes and file mode.

    Never run arbitrary shell code, edit unrelated lines, follow symlinks,
    duplicate the managed block, or overwrite an existing backup.
    """
    validate_bashrc(bashrc)
    existing = bashrc.read_text(encoding="utf-8") if bashrc.exists() else ""
    if BASH_BEGIN in existing:
        return False
    # Respect an older GauMon block, or the manual one-line installation.
    if 'export PATH="$HOME/.local/bin:$PATH"' in existing.splitlines():
        return False
    block = (f"{BASH_BEGIN}\n"
             'case ":$PATH:" in\n'
             '  *":$HOME/.local/bin:"*) ;;\n'
             '  *) export PATH="$HOME/.local/bin:$PATH" ;;\n'
             'esac\n'
             f"{BASH_END}\n")
    # Append-only: never strip existing whitespace or normalize other lines.
    updated = existing + ("\n" if existing and not existing.endswith("\n") else "") + "\n" + block
    mode = stat.S_IMODE(bashrc.stat().st_mode) if bashrc.exists() else 0o644
    if bashrc.exists():
        # Preserve the first original backup forever, use a fresh name thereafter.
        backup = bashrc.with_name(bashrc.name + ".gaumon-backup")
        if backup.exists():
            backup = bashrc.with_name(bashrc.name + f".gaumon-backup-{uuid.uuid4().hex[:12]}")
        shutil.copy2(bashrc, backup, follow_symlinks=False)
    atomic_text(bashrc, updated, mode=mode)
    return True


def copy_application_source(source: Path, stage: Path) -> None:
    """Copy only distribution files; never copy incidental private research data."""
    stage.mkdir(parents=False)
    files = ("run.py", "run.sh", "run_tests.sh", "install.py", "install.sh", "install_user.sh", "patch_source.py",
             "environment.yml", "requirements.txt", "README.md", "LICENSE", "INSTALL_LINUX.md")
    directories = ("gaumon", "gaumon_ext", "assets", "examples", "scripts", "tests")
    for name in files + directories:
        item = source / name
        if not item.exists():
            continue  # validate_app checks the required application core.
        target = stage / name
        if item.is_dir():
            shutil.copytree(item, target, ignore=shutil.ignore_patterns("__pycache__", "*.pyc", ".pytest_cache"))
        else:
            shutil.copy2(item, target)
    # GitHub ZIPs, some mounts, and archives can drop POSIX mode bits.
    # Tests/demo execute mock_g16 directly; always repair its stage copy.
    for name in ("examples/mock_g16", "run.sh", "run_tests.sh"):
        candidate = stage / name
        if candidate.is_file():
            candidate.chmod(candidate.stat().st_mode | 0o111)


def validate_app(app: Path) -> None:
    required = ("run.py", "run.sh", "environment.yml", "requirements.txt", "gaumon",
                "gaumon_ext", "examples", "examples/mock_g16", "tests")
    missing = [name for name in required if not (app / name).exists()]
    if missing:
        fail("Repository checkout incomplete; missing " + ", ".join(missing))
    cli_tests = app / "tests" / "test_cli.py"
    if not cli_tests.is_file():
        fail(f"Missing core CLI regression tests: {cli_tests}")
    cli_test_text = cli_tests.read_text(encoding="utf-8")
    if "test_env_check_reports_expected_name" in cli_test_text:
        fail("Legacy environment-name assertion detected. Apply the current source patch before installing")
    if "test_env_check_reports_runtime_dependencies" not in cli_test_text:
        fail("Unknown CLI environment-check test contract; refusing release installation")
    if '--demo-ui' in (app / 'run.py').read_text(encoding='utf-8'):
        fail('Source still exposes deprecated --demo-ui; apply the v3 patch before installing')
    ui = app / 'gaumon' / 'ui.py'
    if ui.is_file() and 'showFullScreen(' in ui.read_text(encoding='utf-8'):
        fail('Source still uses full-screen GUI; apply the v3 patch before installing')
    for name in ("examples/mock_g16", "run.sh"):
        if not (app / name).is_file() or not os.access(app / name, os.X_OK):
            fail(f"Not executable in staged distribution: {app / name}")


def install(args: argparse.Namespace) -> None:
    if sys.platform != "linux":
        fail("This installation wizard currently supports Linux. No OS changes made.")
    interactive = not args.yes and sys.stdin.isatty()
    conda = select_conda(args.conda_exe)
    mode = pick_mode(args.mode, conda, interactive)
    prefix = args.prefix.expanduser().resolve()
    app = prefix / "app"
    interpreter = environment_python(mode, prefix)
    command = Path.home() / ".local" / "bin" / "gaumon"
    aliases = [command.with_name(name) for name in COMMAND_ALIASES]
    desktop = Path.home() / ".local" / "share" / "applications" / "gaumon.desktop"
    apply_bashrc = args.bashrc or (interactive and not args.no_bashrc and ask_bool("Add ~/.local/bin to ~/.bashrc?", default=True))
    desktop_enabled = not args.no_desktop and (args.desktop or not interactive or ask_bool("Install application-menu shortcut?", default=True))
    log(f"Selected mode: {mode}; environment Python: {interpreter}")
    log(f"Application: {app}; launcher: {command}")
    if app.exists() and not (app.is_dir() and (app / APP_MARKER).is_file()):
        fail(f"Existing application directory is not GauMon-managed: {app}")
    check_managed(command, COMMAND_MARKER)
    for alias in aliases:
        check_managed(alias, COMMAND_MARKER)
    if apply_bashrc:
        validate_bashrc(Path.home() / ".bashrc")
    if desktop_enabled:
        check_managed(desktop, DESKTOP_MARKER)
    if args.dry_run:
        log("DRY RUN: no files changed; no dependencies installed")
        return
    if SOURCE == app or SOURCE in app.parents or SOURCE == prefix:
        fail("Choose an install prefix outside the source checkout")
    prefix.mkdir(parents=True, exist_ok=True)
    stage = prefix / (".app-staging-" + uuid.uuid4().hex[:10])
    backup = prefix / (".app-backup-" + uuid.uuid4().hex[:10])
    moved_old = False
    promoted = False
    try:
        copy_application_source(SOURCE, stage)
        validate_app(stage)
        interpreter = create_environment(mode, prefix, stage, conda, args.python)
        # Make the staged copy behave exactly like the final managed app:
        # run.sh tests must select this new prefix, never an unrelated Conda env.
        (stage / APP_MARKER).write_text(f"mode={mode}\nGauMon managed installation (staging).\n", encoding="utf-8")
        # Resolve Qt imports and genuine bundled unit tests before replacing app.
        subprocess_checked([str(interpreter), "-s", str(stage / "run.py"), "--env-check"], env=clean_python_env())
        subprocess_checked([str(interpreter), "-s", str(stage / "run.py"), "--list-cases"], env=clean_python_env())
        if not args.skip_tests:
            subprocess_checked([str(interpreter), "-s", str(stage / "run.py"), "--test"], env=clean_python_env())
        (stage / APP_MARKER).write_text(f"mode={mode}\nManaged by GauMon install.py.\n", encoding="utf-8")
        if app.exists():
            app.rename(backup)
            moved_old = True
        try:
            stage.rename(app)
            promoted = True
        except Exception:
            if moved_old:
                backup.rename(app)
                moved_old = False
            raise
        atomic_text(command, shell_command(interpreter, app), executable=True)
        for alias in aliases:
            atomic_text(alias, alias_command(command), executable=True)
        if desktop_enabled:
            atomic_text(desktop, desktop_entry(command))
        if apply_bashrc:
            if update_bashrc(Path.home() / ".bashrc"):
                log("Added PATH to ~/.bashrc; backup: ~/.bashrc.gaumon-backup (if previous file existed)")
            else:
                log("~/.bashrc already includes GauMon PATH setup; unchanged")
        if moved_old:
            shutil.rmtree(backup)
        log("INSTALLATION COMPLETE")
        log(f"Verify: {command} --env-check")
        log(f"Launch from any directory: {command}")
        log("Aliases: gaumon, Gaumon, GauMon")
        log(f"PATH setup command (if not already enabled): {bashrc_command()}")
        log("For this existing terminal only: export PATH=\"$HOME/.local/bin:$PATH\" (or open a new shell after opting into ~/.bashrc).")
        log("GUI runtime on a new machine must be verified interactively; Gaussian is NOT installed by GauMon.")
    finally:
        if stage.exists():
            shutil.rmtree(stage)
        if moved_old and not promoted and backup.exists() and not app.exists():
            backup.rename(app)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="GauMon Linux installation wizard (Conda or Python venv)")
    parser.add_argument("--mode", choices=("auto", "conda", "venv"), default="auto", help="environment manager")
    parser.add_argument("--prefix", type=Path, default=Path.home()/".local"/"share"/"gaumon", help="persistent application and environment location")
    parser.add_argument("--conda-exe", help="path to conda executable, if detection fails")
    parser.add_argument("--python", default=sys.executable, help="Python >=3.10 interpreter for venv mode")
    parser.add_argument("--yes", action="store_true", help="noninteractive; select conda if found, otherwise venv")
    bashrc = parser.add_mutually_exclusive_group()
    bashrc.add_argument("--bashrc", action="store_true", help="append a managed PATH block to ~/.bashrc")
    bashrc.add_argument("--no-bashrc", action="store_true", help="never change ~/.bashrc")
    parser.add_argument("--desktop", action="store_true", help="install a GUI desktop shortcut")
    parser.add_argument("--no-desktop", action="store_true", help="do not install a GUI desktop shortcut")
    parser.add_argument("--skip-tests", action="store_true", help="skip unit tests (environment checks always run)")
    parser.add_argument("--dry-run", action="store_true", help="show plan without changing machine")
    args = parser.parse_args(argv)
    try:
        install(args)
        return 0
    except (OSError, RuntimeError, subprocess.CalledProcessError, KeyboardInterrupt) as exc:
        print(f"[GauMon installer] ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
