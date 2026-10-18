import subprocess
import sys
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

class CLITests(unittest.TestCase):
    def call(self,*args):
        return subprocess.run([sys.executable,str(ROOT/'run.py'),*args],cwd=str(ROOT),
                              capture_output=True,text=True,timeout=12)

    def test_deep_stress_cli(self):
        p=self.call('--stress-case','many_steps')
        self.assertEqual(p.returncode,0,p.stderr)
        self.assertIn('8 procedures',p.stdout)
        self.assertIn('complete',p.stdout)

    def test_case_listing(self):
        p=self.call('--list-cases')
        self.assertEqual(p.returncode,0,p.stderr)
        self.assertIn('optfreq',p.stdout)
        self.assertIn('ts',p.stdout)

    def test_good_validation(self):
        p=self.call('--validate',str(ROOT/'examples'/'inputs'/'01_energy_H2.gjf'))
        self.assertEqual(p.returncode,0,p.stderr)

    def test_env_check_reports_expected_name(self):
        p=self.call('--env-check')
        self.assertIn('gaumon-cpu',p.stdout)

    def test_gui_module_compiles_even_if_qt_not_installed(self):
        import py_compile
        py_compile.compile(str(ROOT/'gaumon'/'ui.py'),doraise=True)

    def test_headless_inspect_attached_error(self):
        import json
        fixture=ROOT/'examples'/'error_logs'/'formbx_l103.log'
        result=self.call('--inspect',str(fixture))
        self.assertEqual(result.returncode,0,result.stderr)
        payload=json.loads(result.stdout)
        self.assertEqual(payload['state'],'FAILED')
        self.assertEqual(payload['link'],'L103')
        self.assertIn('FormBX',payload['diagnoses'][0]['title'])

    def test_diagnose_credits_reference(self):
        fixture=ROOT/'examples'/'error_logs'/'formbx_l103.log'
        result=self.call('--diagnose',str(fixture))
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertIn('wongzit.github.io/gaussian-common-errors-and-solutions',result.stdout)
        self.assertIn('19–27–25',result.stdout)
