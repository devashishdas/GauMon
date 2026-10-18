"""Actual Qt smoke tests are skipped if PyQt5/pyqtgraph are unavailable."""
import importlib.util
import os
import unittest
os.environ.setdefault('QT_QPA_PLATFORM','offscreen')
HAS_QT=(importlib.util.find_spec('PyQt5') is not None and
        importlib.util.find_spec('pyqtgraph') is not None)

@unittest.skipUnless(HAS_QT,'PyQt5 and pyqtgraph unavailable in this test runtime')
class QtSmokeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from PyQt5.QtWidgets import QApplication
        cls.app=QApplication.instance() or QApplication([])

    def test_window_without_global_scroll_area(self):
        from gaumon.ui import MainWindow
        from PyQt5.QtWidgets import QScrollArea
        window=MainWindow()
        self.assertEqual(len(window.findChildren(QScrollArea)),0)
        self.assertEqual(window.pages.count(),15)
        self.assertEqual(set(window.plots),{'energy','force','disp','rmsd'})
        self.assertEqual(len(window.eta_tiles),4)
        self.assertEqual(window.link_table.columnCount(),3)
        self.assertFalse(window.nav_buttons[7].isVisible())
        window.open_deep()
        self.assertEqual(window.pages.currentIndex(),7)
        self.assertTrue(window._deep_page_open)
        window.close_deep()
        self.assertFalse(window._deep_page_open)
        self.assertFalse(window.nav_buttons[8].isVisible())
        self.assertFalse(window.nav_buttons[9].isVisible())
        window.open_qcnr()
        self.assertEqual(window.pages.currentIndex(),8)
        self.assertTrue(window._qcnr_page_open)
        self.assertEqual(len(window.qcnr_plots),4)
        window.close_qcnr()
        window.open_freq()
        self.assertEqual(window.pages.currentIndex(),9)
        self.assertTrue(window._freq_page_open)
        self.assertEqual(len(window.freq_charts),2)
        window.close_freq()
        self.assertIsNotNone(window.time_page)
        self.assertIsNotNone(window.input_checker)
        self.assertIsNotNone(window.advanced_eta)
        self.assertIsNotNone(window.warning_page)
        self.assertIsNotNone(window.jobs_page)
        self.assertEqual(window.work_entry.text(),os.getcwd())
        window.show_page(14)
        self.assertEqual(window.pages.currentIndex(),14)
        window.close()

    def test_editor_autoload_and_unsaved_edits(self):
        from gaumon.ui import MainWindow
        from pathlib import Path
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            base=Path(d)
            log=base/'example.log';log.write_text('Entering Gaussian System\n')
            gjf=base/'example.gjf';gjf.write_text('%mem=1GB\n#p opt freq b3lyp/6-31g(d)\n\ntest\n\n0 1\nH 0 0 0\n\n')
            window=MainWindow()
            self.assertTrue(window._autoload_input(log))
            self.assertEqual(window.input_checker.path,gjf.resolve())
            self.assertIn('#p opt freq',window.input_checker.edit.toPlainText())
            self.assertTrue(any(r[1].lower()=='opt' for r in window.input_checker._keyword_rows))
            self.assertTrue(window.input_checker.keyword_table.rowCount()>0)
            window.input_checker.edit.insertPlainText(' edited')
            original=window.input_checker.edit.toPlainText()
            self.assertTrue(window._autoload_input(log))
            self.assertEqual(window.input_checker.edit.toPlainText(),original)
            window.close()


    def test_saved_session_records_do_not_crash_first_draw(self):
        from gaumon.ui import MainWindow
        from gaumon_ext.session import SessionStore
        from pathlib import Path
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            window=MainWindow()
            try:
                log=Path(d)/'saved_without_reader.log'
                log.write_text('Entering Gaussian System\n')
                window.hub.store=SessionStore(Path(d)/'session.json')
                window.hub.store.add(log)
                # A saved job has no active parser snapshot before explicit resume.
                window.jobs_page.redraw()
                self.assertEqual(window.jobs_page.table.rowCount(),1)
                self.assertEqual(window.jobs_page.table.item(0,1).text(),
                                 'SAVED · RESUME TO ATTACH')
                self.assertEqual(window.jobs_page.table.item(0,3).text(),'—')
                window.hub.resume(window.hub.store.items[0]['id'])
                window.jobs_page.redraw()
                self.assertEqual(window.jobs_page.table.rowCount(),1)
                self.assertIn('UNVERIFIED',window.jobs_page.table.item(0,1).text())
            finally:
                window.close()

    def test_testcase_panel_has_four_rows(self):
        from gaumon.ui import MainWindow
        window=MainWindow()
        self.assertEqual(window.test_table.rowCount(),7)
        window.close()
