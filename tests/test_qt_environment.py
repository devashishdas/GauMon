"""Dependency diagnostics must fail when a broken Qt installation is on sys.path."""
import unittest
from unittest.mock import patch
from gaumon_ext.qtcheck import qt_environment_report, REQUIRED


class QtEnvironmentTests(unittest.TestCase):
    def test_bad_sip_is_not_reported_ready(self):
        def importer(name):
            if name == 'PyQt5.sip':
                raise ModuleNotFoundError("No module named 'PyQt5.sip'")
            if name == 'PyQt5.QtCore':
                class Core:
                    QT_VERSION_STR = '5.15.0'
                    PYQT_VERSION_STR = '5.15.0'
                return Core
            return type('Mod', (), {'__file__': f'/test/{name}.so'})
        with patch('gaumon_ext.qtcheck.importlib.import_module', side_effect=importer):
            good, report = qt_environment_report(verbose=False)
        self.assertFalse(good)
        self.assertIn("FAILED  PyQt5.sip", report)
        self.assertIn('pyqt5-sip', report)

    def test_all_imports_check_the_actual_modules(self):
        seen=[]
        def importer(name):
            seen.append(name)
            if name == 'PyQt5.QtCore':
                class Core:
                    QT_VERSION_STR = '5.15.0'
                    PYQT_VERSION_STR = '5.15.0'
                return Core
            return type('Mod', (), {'__file__': f'/test/{name}.so'})
        with patch('gaumon_ext.qtcheck.importlib.import_module', side_effect=importer):
            good, report = qt_environment_report(verbose=False)
        self.assertTrue(good)
        self.assertEqual(seen, list(REQUIRED))
        self.assertIn('PyQt5.QtWidgets',report)

if __name__ == '__main__':
    unittest.main()
