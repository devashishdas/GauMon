"""Security tests for the consolidated v7 desktop, not the legacy GauMon.py."""
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class SafetyTests(unittest.TestCase):
    def test_no_global_kill(self):
        sources=[ROOT/'gaumon'/'engine.py',ROOT/'gaumon'/'ui.py',ROOT/'gaumon_ext'/'gui.py']
        for src in sources:
            with self.subTest(src=src.name):
                text=src.read_text()
                self.assertNotIn('pkill -9',text)
                self.assertNotIn('os.system(',text)
        text=sources[0].read_text()
        self.assertIn("if self.mode != 'managed'",text)
        self.assertIn('os.killpg(self.pgid',text)
    def test_launcher_preserves_original_working_directory(self):
        content=(ROOT/'run.sh').read_text()
        self.assertNotIn('\ncd "$(dirname "${BASH_SOURCE[0]}")"',content)
        self.assertIn('python "$ROOT/run.py"',content)
    def test_desktop_does_not_autoresume(self):
        text=(ROOT/'run.py').read_text()
        self.assertNotIn('hub.resume(',text)
        self.assertNotIn('monitor.launch(', (ROOT/'gaumon_ext'/'session.py').read_text())
