"""Static release contracts: no demo command, no immersive full-screen startup."""
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]

class ReleaseBehaviorTests(unittest.TestCase):
    def test_demo_ui_removed_from_cli(self):
        source = (ROOT / 'run.py').read_text()
        self.assertNotIn('--demo-ui', source)
        self.assertNotIn('a.demo_ui', source)

    def test_gui_is_maximized_and_not_fullscreen(self):
        source = (ROOT / 'gaumon/ui.py').read_text()
        self.assertIn('window.showMaximized()', source)
        self.assertNotIn('window.showFullScreen()', source)
        self.assertNotIn('toggle_fullscreen', source)
        self.assertNotIn('self.full_btn', source)
        self.assertNotIn('mock_on_start', source)

if __name__ == '__main__':
    unittest.main()
