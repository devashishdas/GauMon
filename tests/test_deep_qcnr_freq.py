"""Distinguish QCNR linear solves from vibrational frequency modes."""
import math
from pathlib import Path
import subprocess
import sys
import unittest
from gaumon.engine import Parser

ROOT=Path(__file__).resolve().parents[1]
STRESS=ROOT/'examples'/'stress_logs'

class NestedSolverTests(unittest.TestCase):
    def test_complete_and_interrupted_solve(self):
        p=Parser();p.feed((STRESS/'qcnr_lineq1.log').read_text())
        meta=p.snapshot()['deep_qcnr']
        self.assertEqual(len(meta['solves']),2)
        self.assertEqual(meta['outer_iterations'],2)
        first=p.qcnr.get_solve(1)
        second=p.qcnr.get_solve(2)
        self.assertEqual(first['status'],'converged')
        self.assertEqual(first['actual_convergence']['iterations_reported'],3)
        self.assertEqual(first['rows'][-1]['noncon'],0)
        self.assertAlmostEqual(first['rows'][-1]['rms'],1.82e-9)
        self.assertEqual(second['status'],'running')
        self.assertEqual(second['rows'][0]['iteration'],0)
        self.assertEqual(second['rows'][-1]['iteration'],2)
        self.assertEqual(p.qcnr.outer_history()[1]['delta_e'],-0.001432882648)
        self.assertEqual(len(p.snapshot()['deep_freq']['groups']),0)

    def test_split_file_chunks_preserve_history(self):
        text=(STRESS/'qcnr_lineq1.log').read_text()
        a=Parser();a.feed(text)
        b=Parser()
        for i in range(0,len(text),3):b.feed(text[i:i+3])
        self.assertEqual(a.snapshot()['deep_qcnr'],b.snapshot()['deep_qcnr'])
        self.assertEqual(a.qcnr.get_solve(1),b.qcnr.get_solve(1))

    def test_restarted_solver_at_same_outer_iteration(self):
        text=''' (Enter /g/l508.exe)
 Quadratic Convergence SCF Method.
 Iteration 1 EE= -200.0 Grad=2.00D-03 Max rot=0.5 deg.
 LinEq1: Iter=0 NonCon=1 RMS=1.0D-02 Max=1.0D-01 NDo=1
 LinEq1: Iter=1 NonCon=1 RMS=5.0D-03 Max=5.0D-02 NDo=1
 LinEq1: Iter=0 NonCon=1 RMS=2.0D-02 Max=2.0D-01 NDo=1
 LinEq1: Iter=1 NonCon=0 RMS=1.0D-09 Max=1.0D-08 NDo=1
'''
        p=Parser();p.feed(text)
        solves=p.qcnr.overview()['solves']
        self.assertEqual([a['status'] for a in solves],['interrupted','converged'])
        self.assertEqual([a['iterations'] for a in solves],[2,2])

    def test_end_of_log_unterminated_line_is_not_fabricated(self):
        p=Parser()
        p.feed(' (Enter /g/l508.exe)\n Quadratic Convergence SCF Method.\n'
               ' Iteration 1 EE= -200 Grad=1D-03 Max rot=0.1 deg.\n'
               ' LinEq1: Iter=0 NonCon=1 RMS=1D-02 Max=1D-01 NDo=1\n'
               ' LinEq1: Iter=1 NonCon=0 RMS=')
        self.assertEqual(p.qcnr.overview()['solves'][0]['iterations'],1)
        self.assertEqual(p.qcnr.get_solve()['status'],'running')

    def test_converged_inner_solver_does_not_prove_scf_converged(self):
        p=Parser()
        p.feed((STRESS/'qcnr_lineq1.log').read_text())
        self.assertIsNone(p.snapshot()['termination'])
        self.assertEqual(p.qcnr.get_solve(1)['status'],'converged')
        self.assertFalse(p.snapshot()['opt_completed'])

class FrequencyModeTests(unittest.TestCase):
    def test_real_frequency_groups(self):
        p=Parser();p.feed((STRESS/'harmonic_modes.log').read_text())
        meta=p.snapshot()['deep_freq']
        self.assertEqual(len(meta['groups']),1)
        group=p.freq_deep.get_group()
        self.assertEqual(group['status'],'printed')
        self.assertEqual(len(group['modes']),6)
        self.assertEqual([r['frequency'] for r in group['modes'][:3]],[-120.52,230.12,456.78])
        self.assertEqual(group['modes'][0]['red_mass'],1.23)
        self.assertEqual(group['modes'][0]['force_constant'],0.01)
        self.assertEqual(group['modes'][0]['ir_intensity'],23.1)
        self.assertEqual(group['modes'][-1]['ir_intensity'],3.0)
        self.assertEqual(meta['groups'][0]['imaginary'],1)
        self.assertEqual(len(p.snapshot()['deep_qcnr']['solves']),0)

    def test_frequencies_chunks_and_separate_group(self):
        text=(STRESS/'harmonic_modes.log').read_text()
        a=Parser();a.feed(text)
        b=Parser()
        for i in range(0,len(text),2):b.feed(text[i:i+2])
        self.assertEqual(a.snapshot()['deep_freq'],b.snapshot()['deep_freq'])
        self.assertEqual(a.freq_deep.get_group(),b.freq_deep.get_group())
        b.feed(' Harmonic frequencies (cm**-1)\n Frequencies -- 50.0\n')
        self.assertEqual(len(b.snapshot()['deep_freq']['groups']),2)
        self.assertEqual(b.freq_deep.get_group()['modes'][0]['frequency'],50.0)

    def test_cli_reports(self):
        for switch, filename, expected in (('--qcnr-report','qcnr_lineq1.log','Linear equation'),
                                            ('--freq-report','harmonic_modes.log','6 modes')):
            a=subprocess.run([sys.executable,str(ROOT/'run.py'),switch,str(STRESS/filename)],
                              cwd=ROOT, capture_output=True,text=True,timeout=12)
            self.assertEqual(a.returncode,0,a.stderr)
            self.assertIn(expected.lower(),a.stdout.lower())

if __name__=='__main__':unittest.main()
