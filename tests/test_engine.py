import os
import tempfile
import time
import unittest
from pathlib import Path
import numpy as np
from gaumon.engine import Parser,Monitor,kabsch
from gaumon.inputs import write_case

ROOT=Path(__file__).resolve().parents[1]
MOCK=ROOT/'examples'/'mock_g16'

class TestParser(unittest.TestCase):
    def test_kabsch_rigid_alignment(self):
        coords=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
        moved=coords@np.array([[0.,1.,0.],[-1.,0.,0.],[0.,0.,1.]])+[5.,3.,8.]
        rmsd,maxdisp=kabsch(coords,moved)
        self.assertLess(rmsd,1e-10)
        self.assertLess(maxdisp,1e-10)

    def test_kabsch_rejects_mismatched_count(self):
        self.assertEqual(kabsch([[0,0,0]],[[0,0,0],[1,0,0]]),(None,None))

    def test_incremental_partial_lines(self):
        parser=Parser()
        parser.feed('SCF Done:  E(RB3LYP) =  -10.000000 A.U. after  6')
        self.assertIsNone(parser.snapshot()['energy'])
        parser.feed(' cycles\n')
        self.assertEqual(parser.snapshot()['energy'],-10.)
        self.assertEqual(parser.snapshot()['cycles'],6)

    def test_frequency_and_optimization(self):
        parser=Parser()
        parser.feed('Step number 1 out of a maximum of 40\nOptimization completed.\n'
                    'Frequencies --  -850.2 230.0 500.0\n')
        snap=parser.snapshot()
        self.assertEqual(snap['step'],1)
        self.assertTrue(snap['opt_completed'])
        self.assertEqual(snap['imaginary_modes'],1)

    def test_rescan_for_link_state(self):
        parser=Parser();parser.feed('Error termination\nEntering Gaussian System\n')
        self.assertIsNone(parser.termination)


class TestMockCases(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.work=Path(self.tmp.name)
        self.monitor=Monitor()

    def run_mock(self,key,delay='0.02'):
        inp=write_case(key,self.work)
        old=os.environ.get('GAUMON_MOCK_STEP_SECONDS')
        os.environ['GAUMON_MOCK_STEP_SECONDS']=delay
        try:
            self.monitor.launch(str(MOCK),str(inp),str(self.work),key+'.log')
            deadline=time.monotonic()+12
            while time.monotonic()<deadline:
                s=self.monitor.poll()
                if s['state'] in self.monitor.TERMINAL:
                    return s
                time.sleep(0.08)
            self.fail('mock subprocess did not reach terminal status')
        finally:
            if old is None:os.environ.pop('GAUMON_MOCK_STEP_SECONDS',None)
            else:os.environ['GAUMON_MOCK_STEP_SECONDS']=old
            if self.monitor.proc and self.monitor.proc.poll() is None:
                self.monitor.stop(force=True)

    def test_single_point_energy(self):
        s=self.run_mock('sp')
        self.assertEqual(s['state'],'COMPLETED')
        self.assertEqual(len(s['scf']),1)
        self.assertEqual(s['frequencies'],0)

    def test_optimization_plus_frequency(self):
        s=self.run_mock('optfreq')
        self.assertEqual(s['state'],'COMPLETED')
        self.assertTrue(s['opt_completed'])
        self.assertGreaterEqual(s['frequencies'],1)
        self.assertGreaterEqual(s['geometries'],1)
        self.assertEqual(s['imaginary_modes'],0)

    def test_frequency_only(self):
        s=self.run_mock('freq')
        self.assertEqual(s['state'],'COMPLETED')
        self.assertFalse(s['opt_completed'])
        self.assertGreaterEqual(s['frequencies'],1)

    def test_transition_state_mock_has_one_negative_frequency(self):
        s=self.run_mock('ts')
        self.assertEqual(s['state'],'COMPLETED')
        self.assertTrue(s['opt_completed'])
        self.assertEqual(s['imaginary_modes'],1)
        self.assertLess(min(s['frequency_values']),0)

    def test_timer_freezes_after_termination(self):
        s=self.run_mock('sp');elapsed=s['elapsed']
        time.sleep(.07)
        self.assertEqual(self.monitor.poll()['elapsed'],elapsed)

    def test_stop_only_owned_mock(self):
        inp=write_case('ts',self.work)
        os.environ['GAUMON_MOCK_STEP_SECONDS']='3.0'
        try:
            self.monitor.launch(str(MOCK),str(inp),str(self.work),'stop_test.log')
            self.monitor.stop(force=True)
            deadline=time.monotonic()+8
            while time.monotonic()<deadline:
                s=self.monitor.poll()
                if s['state']=='STOPPED':break
                time.sleep(.05)
            self.assertEqual(s['state'],'STOPPED')
        finally:
            os.environ.pop('GAUMON_MOCK_STEP_SECONDS',None)
            if self.monitor.proc and self.monitor.proc.poll() is None:
                self.monitor.stop(force=True)

    def test_attached_completed_log_freezes(self):
        log=self.work/'old.log'
        log.write_text(' SCF Done:  E(RHF) =  -1.00     A.U. after  5 cycles\nNormal termination of Gaussian 16\n')
        s=self.monitor.attach(str(log))
        done=self.monitor.poll()
        self.assertEqual(done['state'],'COMPLETED')
        elapsed=done['elapsed']
        time.sleep(.06)
        self.assertEqual(self.monitor.poll()['elapsed'],elapsed)

    def test_external_process_cannot_be_stopped(self):
        log=self.work/'ext.log';log.write_text('Starting\n')
        self.monitor.attach(str(log))
        with self.assertRaises(ValueError):self.monitor.stop(force=True)


    def test_eta_and_link_history_are_generated(self):
        s=self.run_mock('optfreq')
        self.assertTrue(s['eta_records'])
        self.assertEqual([r['extra_cycles'] for r in s['eta_records']], [10,50,100,150])
        self.assertGreaterEqual(len(s['step_history']), 1)
        self.assertGreaterEqual(len(s['link_history']), 1)
        self.assertIn('description', s['link_history'][0])

    def test_refuse_nonempty_log_overwrite(self):
        inp=write_case('sp',self.work)
        (self.work/'sp.log').write_text('valuable old output')
        with self.assertRaises(ValueError):
            self.monitor.launch(str(MOCK),str(inp),str(self.work),'sp.log')
        self.assertEqual((self.work/'sp.log').read_text(),'valuable old output')
