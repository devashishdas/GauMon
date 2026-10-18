"""Tests for consolidated v7 execution pathways, no Gaussian or Qt needed."""
from __future__ import annotations
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from gaumon_ext.analytics import IncrementalLog
from gaumon_ext.checker import check_input
from gaumon_ext.session import SessionStore

ROOT=Path(__file__).resolve().parents[1]

class ConsolidationTests(unittest.TestCase):
    def test_current_directory_kept_for_relative_log(self):
        with tempfile.TemporaryDirectory() as d:
            file=Path(d)/'sample.log'
            file.write_text('Leave Link 502 at Sun Oct 18 11:00:00 2026, MaxMem= 123 cpu: 20.2 elap: 5.1\n')
            env=os.environ.copy();env['CONDA_DEFAULT_ENV']='gaumon-cpu'
            proc=subprocess.run(['bash',str(ROOT/'run.sh'),'--timing-report','sample.log'],
                                cwd=d,env=env,capture_output=True,text=True,timeout=15)
            self.assertEqual(proc.returncode,0,proc.stderr)
            self.assertIn('5.10',proc.stdout)
            self.assertIn('L502',proc.stdout)
    def test_timing_report_against_synthetic_fixture(self):
        run=subprocess.run([sys.executable,str(ROOT/'run.py'),'--timing-report',
                            str(ROOT/'examples/mock_timing.log')],capture_output=True,text=True,timeout=12)
        self.assertEqual(run.returncode,0,run.stderr)
        self.assertIn('4 completed Link calls',run.stdout)
        self.assertIn('17.00',run.stdout)
    def test_checker_cli_blocks_bad_structure(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'bad.gjf';p.write_text('#p hf/3-21g sp\nNO BLANK TITLE\n0 1\nH 0 0 0\n')
            run=subprocess.run([sys.executable,str(ROOT/'run.py'),'--check-input',str(p)],
                                capture_output=True,text=True,timeout=12)
            self.assertEqual(run.returncode,1)
            self.assertIn('ERROR',run.stdout)
    def test_saved_job_requires_explicit_attachment(self):
        with tempfile.TemporaryDirectory() as folder:
            d=Path(folder)
            log=d/'job.log';log.write_text('Gaussian job log\n')
            store=SessionStore(d/'session.json');store.add(log)
            restored=SessionStore(d/'session.json')
            self.assertEqual(len(restored.items),1)
            self.assertIsNone(restored.items[0]['pid'])
            self.assertEqual(restored.items[0]['last_known'],'UNVERIFIED')
            self.assertIsNone(IncrementalLog(log).snapshot.last_growth)
    def test_advanced_modules_all_exist(self):
        for name in ('deep_scf.py','deep_qcnr.py','deep_freq.py'):
            self.assertTrue((ROOT/'gaumon'/name).is_file())
        self.assertIn('TimePage(self.hub)',(ROOT/'gaumon/ui.py').read_text())
        self.assertIn('JobsPage(self.hub',(ROOT/'gaumon/ui.py').read_text())

if __name__ == '__main__':unittest.main()
