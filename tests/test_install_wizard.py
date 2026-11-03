"""Installer tests use temporary homes; no conda, pip, GUI, or Gaussian required."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest import mock

INSTALLER = Path(__file__).resolve().parents[1] / "install.py"
spec = importlib.util.spec_from_file_location("gaumon_installer", INSTALLER)
installer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(installer)


class WizardTests(unittest.TestCase):
    def test_managed_bashrc_idempotent(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / ".bashrc"
            path.write_text("# existing config\n", encoding="utf-8")
            self.assertTrue(installer.update_bashrc(path))
            self.assertFalse(installer.update_bashrc(path))
            content = path.read_text()
            self.assertEqual(content.count(installer.BASH_BEGIN), 1)
            self.assertIn('# existing config', content)
            self.assertEqual((Path(d) / ".bashrc.gaumon-backup").read_text(), "# existing config\n")

    def test_incomplete_block_refuses_edit(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / ".bashrc"
            path.write_text(installer.BASH_BEGIN + "\n")
            with self.assertRaises(RuntimeError):
                installer.update_bashrc(path)

    def test_launcher_preserves_arguments_without_shell_expansion(self):
        python = Path('/tmp/env with space/bin/python')
        app = Path('/tmp/repo with space')
        text = installer.shell_command(python, app)
        self.assertIn(installer.COMMAND_MARKER, text)
        self.assertIn('"$@"', text)
        self.assertIn("'/tmp/env with space/bin/python'", text)

    def test_legacy_gaumon_command_is_upgradable(self):
        with tempfile.TemporaryDirectory() as d:
            launcher = Path(d) / "gaumon"
            launcher.write_text("#!/usr/bin/env bash\n# GAUMON_USER_LAUNCHER\n")
            installer.check_managed(launcher, installer.COMMAND_MARKER)

    def test_existing_manual_path_is_not_duplicated(self):
        with tempfile.TemporaryDirectory() as d:
            bashrc = Path(d) / ".bashrc"
            line = 'export PATH="$HOME/.local/bin:$PATH"'
            bashrc.write_text(line + "\n")
            self.assertFalse(installer.update_bashrc(bashrc))
            self.assertEqual(bashrc.read_text().count(line), 1)

    def test_unmanaged_command_refused(self):
        with tempfile.TemporaryDirectory() as d:
            file = Path(d) / 'gaumon'
            file.write_text('unrelated executable')
            with self.assertRaises(RuntimeError):
                installer.check_managed(file, installer.COMMAND_MARKER)

    def test_detect_conda_explicit_missing(self):
        self.assertIsNone(installer.select_conda('/not/a/real/conda'))

    def test_conda_requires_conda(self):
        with self.assertRaises(RuntimeError):
            installer.pick_mode('conda', None, interactive=False)

    def test_dry_run_writes_nothing(self):
        with tempfile.TemporaryDirectory() as d:
            home = Path(d) / 'home'
            home.mkdir()
            prefix = home / 'app prefix'
            with mock.patch.object(installer.Path, 'home', return_value=home):
                with mock.patch('sys.stdin.isatty', return_value=False):
                    code = installer.main(['--mode', 'venv', '--prefix', str(prefix), '--yes', '--dry-run', '--no-desktop'])
            self.assertEqual(code, 0)
            self.assertFalse(prefix.exists())

    def test_end_to_end_install_with_mock_environment(self):
        # Full copy/promote/launcher/bashrc/desktop lifecycle without pip/network.
        import os
        import subprocess
        import sys
        with tempfile.TemporaryDirectory() as d:
            base = Path(d)
            home = base / "home"
            home.mkdir()
            source = base / "source"
            source.mkdir()
            for file in ("environment.yml", "requirements.txt"):
                (source / file).write_text("# synthetic test file\n")
            for folder in ("gaumon", "gaumon_ext", "examples", "tests"):
                (source / folder).mkdir()
            (source / "run.py").write_text(
                "import sys\nprint('MOCK_GAUMON', *sys.argv[1:])\n"
            )
            (source / "run.sh").write_text('#!/usr/bin/env bash\necho MOCK_RUNSH\n')
            (source / "tests/test_cli.py").write_text(
                "    def test_env_check_reports_runtime_dependencies(self):\n"
                "        pass\n"
            )
            (source / "examples/mock_g16").write_text('#!/usr/bin/env python3\nprint("SYNTHETIC_G16")\n')
            # Deliberately remove executable bits. The installer must restore
            # them before it tests or promotes the staged application.
            (source / "run.sh").chmod(0o644)
            (source / "examples/mock_g16").chmod(0o644)
            prefix = home / "gaumon-location"
            with mock.patch.object(installer.Path, "home", return_value=home), \
                 mock.patch.object(installer, "SOURCE", source), \
                 mock.patch.object(installer, "create_environment", return_value=Path(sys.executable)):
                status = installer.main([
                    "--mode", "venv", "--prefix", str(prefix),
                    "--yes", "--bashrc", "--desktop",
                ])
                self.assertEqual(status, 0)
                launcher = home / ".local/bin/gaumon"
                self.assertTrue(launcher.exists())
                result = subprocess.run([str(launcher), "--version"], capture_output=True, text=True,
                                        env={**os.environ, "PYTHONPATH": "/tmp/evil"})
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("MOCK_GAUMON --version", result.stdout)
                for alias_name in installer.COMMAND_ALIASES:
                    alias = launcher.with_name(alias_name)
                    self.assertTrue(alias.is_file(), alias_name)
                    called = subprocess.run([str(alias), "--list-cases"], cwd=str(base),
                                            capture_output=True, text=True)
                    self.assertEqual(called.returncode, 0, called.stderr)
                    self.assertIn("MOCK_GAUMON --list-cases", called.stdout)
                self.assertTrue((prefix / "app" / installer.APP_MARKER).is_file())
                self.assertTrue((prefix / "app" / "run.sh").is_file())
                self.assertTrue((prefix / "app" / "run.sh").stat().st_mode & 0o111)
                self.assertTrue((prefix / "app" / "examples/mock_g16").stat().st_mode & 0o111)
                demo = subprocess.run([str(prefix / "app" / "examples/mock_g16")],
                                      capture_output=True, text=True)
                self.assertEqual(demo.returncode, 0, demo.stderr)
                self.assertIn("SYNTHETIC_G16", demo.stdout)
                self.assertIn(installer.BASH_BEGIN, (home / ".bashrc").read_text())
                self.assertIn(installer.DESKTOP_MARKER, (home / ".local/share/applications/gaumon.desktop").read_text())

    def test_validate_app_requires_prefix_safe_environment_test(self):
        import os
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            for name in ("run.py", "run.sh", "environment.yml", "requirements.txt"):
                (root / name).write_text("# fixture\n")
            for name in ("gaumon", "gaumon_ext", "examples", "tests"):
                (root / name).mkdir()
            (root / "gaumon/ui.py").write_text("# window maximized\n")
            (root / "examples/mock_g16").write_text("#!/usr/bin/env python3\n")
            (root / "run.sh").chmod(0o755)
            (root / "examples/mock_g16").chmod(0o755)
            (root / "tests/test_cli.py").write_text(
                "    def test_env_check_reports_expected_name(self):\n"
                "        self.assertIn('gaumon-cpu', 'some output')\n"
            )
            with self.assertRaisesRegex(RuntimeError, "Legacy environment-name assertion"):
                installer.validate_app(root)
            (root / "tests/test_cli.py").write_text(
                "    def test_env_check_reports_runtime_dependencies(self):\n"
                "        pass\n"
            )
            installer.validate_app(root)

    def test_staging_copies_run_sh_and_repairs_mock_permissions(self):
        import os
        import subprocess
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root / "examples").mkdir()
            (root / "run.sh").write_text("#!/usr/bin/env bash\necho started\n")
            (root / "examples/mock_g16").write_text("#!/usr/bin/env python3\nprint('mock ok')\n")
            for name in ("run.sh", "examples/mock_g16"):
                (root / name).chmod(0o644)
            stage = root / "stage"
            installer.copy_application_source(root, stage)
            self.assertTrue(os.access(stage / "run.sh", os.X_OK))
            self.assertTrue(os.access(stage / "examples/mock_g16", os.X_OK))
            output = subprocess.run([str(stage / "examples/mock_g16")],
                                    capture_output=True, text=True)
            self.assertEqual(output.returncode, 0, output.stderr)
            self.assertIn("mock ok", output.stdout)

    def test_validate_app_refuses_missing_run_sh(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            for name in ("run.py", "environment.yml", "requirements.txt"):
                (root / name).touch()
            for name in ("gaumon", "gaumon_ext", "examples", "tests"):
                (root / name).mkdir()
            (root / "examples/mock_g16").touch(mode=0o755)
            with self.assertRaisesRegex(RuntimeError, "run.sh"):
                installer.validate_app(root)

    def test_clean_environment_blocks_qt_and_python_conflicts(self):
        with mock.patch.dict('os.environ', {'PYTHONPATH': '/old', 'CONDA_PREFIX': '/old/env',
                                          'QT_QPA_PLATFORM_PLUGIN_PATH': '/old/qt'}):
            child_env = installer.clean_python_env()
        for key in ('PYTHONPATH', 'CONDA_PREFIX', 'QT_QPA_PLATFORM_PLUGIN_PATH'):
            self.assertNotIn(key, child_env)
        self.assertEqual(child_env['PYTHONNOUSERSITE'], '1')

    def test_managed_run_sh_prefers_selected_mode(self):
        import os
        import subprocess
        with tempfile.TemporaryDirectory() as d:
            prefix = Path(d)
            app = prefix / "app"
            app.mkdir()
            (app / "run.sh").write_text((INSTALLER.parent / "run.sh").read_text())
            (app / ".gaumon-managed-app").write_text("mode=venv\n")
            (app / "run.py").write_text("print('ignore')\n")
            (prefix / "venv/bin").mkdir(parents=True)
            (prefix / "conda-env/bin").mkdir(parents=True)
            for mode, label in (("venv", "CORRECT_VENV"), ("conda-env", "WRONG_CONDA")):
                script = prefix / mode / "bin/python"
                script.write_text(f"#!/usr/bin/env bash\necho {label}\n")
                script.chmod(0o755)
            result = subprocess.run(["bash", str(app / "run.sh"), "--version"],
                                    env={**os.environ, "CONDA_PREFIX": "/wrong/active/conda"},
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("CORRECT_VENV", result.stdout)
            self.assertNotIn("WRONG_CONDA", result.stdout)

    def test_desktop_entry_generated(self):
        text = installer.desktop_entry(Path('/home/user/.local/bin/gaumon'))
        self.assertIn('Exec="/home/user/.local/bin/gaumon" --gui', text)
        self.assertIn(installer.DESKTOP_MARKER, text)

    def test_environment_interpreter_paths(self):
        root = Path('/tmp/gaumon')
        self.assertEqual(installer.environment_python('conda', root), root / 'conda-env/bin/python')
        self.assertEqual(installer.environment_python('venv', root), root / 'venv/bin/python')


if __name__ == '__main__':
    unittest.main()

# Added regression checks for PATH integrity and mixed-case entry points.
class PathSafetyV3Tests(unittest.TestCase):
    def test_bashrc_symlink_is_never_modified(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            original=root/'original'
            original.write_text('unaltered\n')
            link=root/'.bashrc'
            link.symlink_to(original)
            with self.assertRaisesRegex(RuntimeError, 'symlinked'):
                installer.update_bashrc(link)
            self.assertEqual(original.read_text(), 'unaltered\n')

    def test_bashrc_preserves_original_permissions_and_content(self):
        import stat
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'.bashrc'
            original='if test -f /tmp/a; then\n  echo run\nfi\n'
            path.write_text(original)
            path.chmod(0o600)
            self.assertTrue(installer.update_bashrc(path))
            self.assertTrue(path.read_text().startswith(original))
            self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o600)
            self.assertIn('case ":$PATH:" in',path.read_text())
            self.assertFalse(installer.update_bashrc(path))

    def test_preexisting_backup_not_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'.bashrc'
            path.write_text('original\n')
            backup=Path(d)/'.bashrc.gaumon-backup'
            backup.write_text('previous-backup')
            installer.update_bashrc(path)
            self.assertEqual(backup.read_text(), 'previous-backup')
            self.assertEqual(len(list(Path(d).glob('.bashrc.gaumon-backup-*'))), 1)

    def test_alias_exec_shell_quoting(self):
        script=installer.alias_command(Path('/tmp/a path/gaumon'))
        self.assertIn(installer.COMMAND_MARKER, script)
        self.assertIn("exec '/tmp/a path/gaumon' \"$@\"", script)
        self.assertEqual(installer.COMMAND_ALIASES, ('Gaumon','GauMon'))
