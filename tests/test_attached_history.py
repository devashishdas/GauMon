"""Regression tests for attached log replay, update stamps, ETA confidence, and FormBX."""
import tempfile
import time
import unittest
from pathlib import Path
from gaumon.diagnostics import diagnose_log, GUIDE_URL
from gaumon.engine import Monitor

ROOT=Path(__file__).resolve().parents[1]


class AttachedRegressionTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.log=Path(self.tmp.name)/'example.log'

    def test_formbx_specific_before_generic_and_source_credit(self):
        text=(ROOT/'examples'/'error_logs'/'formbx_l103.log').read_text()
        diagnoses=diagnose_log(text)
        self.assertIn('FormBX', diagnoses[0].title)
        self.assertIn('19–27–25', diagnoses[0].description)
        self.assertIn('20–19–27–25', diagnoses[0].description)
        self.assertIn('Opt=Cartesian', diagnoses[0].action)
        self.assertEqual(diagnoses[0].source_url,GUIDE_URL)
        self.log.write_text(text)
        m=Monitor();m.attach(self.log)
        s=m.poll()
        self.assertEqual(s['state'],'FAILED')
        self.assertIn('FormBX', s['diagnostic_text'])
        self.assertEqual(s['link'],'L103')

    def test_attach_replays_every_link_and_step_without_false_timing(self):
        self.log.write_text(''' (Enter /software/g16/l103.exe)
 Step number 1 out of a maximum of 108
 SCF Done:  E(RB3LYP) =  -10.000 A.U. after 14 cycles
 (Enter /software/g16/l502.exe)
 (Enter /software/g16/l103.exe)
 Step number 2 out of a maximum of 108
 SCF Done:  E(RB3LYP) =  -11.000 A.U. after 11 cycles
 (Enter /software/g16/l502.exe)
 (Enter /software/g16/l103.exe)
 Step number 3 out of a maximum of 108
''')
        m=Monitor();m.attach(self.log)
        s=m.poll()
        self.assertEqual(s['state'],'ATTACHED')  # log-only, process unverified
        self.assertEqual(len(s['step_history']),3)
        self.assertEqual(len(s['link_history']),5)
        self.assertEqual([r['link'] for r in s['link_history']],
                         ['L103','L502','L103','L502','L103'])
        self.assertEqual(s['cycle_sample_count'],0)
        self.assertIsNone(s['eta_records'][0]['seconds'])
        self.assertEqual(s['scf'][-1][2],11)
        self.assertEqual(s['link_description'],'Berny optimizations to minima and TS, STQN transition state searches')
        self.assertTrue(s['revision'] > 0)
        self.assertTrue(all(e['cycle_time_s'] is None for e in s['step_history']))

    def test_attach_timestamp_from_gaussian_is_used_when_available(self):
        self.log.write_text(''' Leave Link 103 at Sun Oct 18 04:01:00 2026, cpu: 1.0 elap: 1.0
 Step number 1 out of a maximum of 30
 Leave Link 103 at Sun Oct 18 04:01:20 2026, cpu: 1.0 elap: 1.0
 Step number 2 out of a maximum of 30
 Leave Link 103 at Sun Oct 18 04:01:45 2026, cpu: 1.0 elap: 1.0
 Step number 3 out of a maximum of 30
''')
        m=Monitor();m.attach(self.log)
        s=m.poll()
        self.assertEqual(s['cycle_sample_count'],2)
        self.assertEqual(s['cycle_durations'],[20.0,25.0])
        self.assertEqual(s['cycle_timing_source'],'Gaussian link-time proxy')
        self.assertEqual(s['eta_records'][0]['seconds'],225.0)

    def test_attached_tail_refreshes_after_130_lines(self):
        self.log.write_text(''.join(f'line {i}\n' for i in range(220)))
        m=Monitor();m.attach(self.log)
        first=m.poll()
        self.assertEqual(len(first['tail'].splitlines()),130)
        old_revision=first['revision']
        with self.log.open('a') as f:f.write('further log output appeared\n')
        time.sleep(.38)
        after=m.poll()
        self.assertGreater(after['revision'],old_revision)
        self.assertEqual(len(after['tail'].splitlines()),130)
        self.assertIn('further log output appeared',after['tail'])

    def test_attach_new_live_cycles_without_history_times(self):
        self.log.write_text('Step number 1 out of a maximum of 40\n')
        m=Monitor();m.attach(self.log)
        m.poll()
        time.sleep(.4)
        with self.log.open('a') as f:f.write('Step number 2 out of a maximum of 40\n')
        m.poll()
        # First post-attachment cycle is not trusted as a complete interval.
        self.assertEqual(m.snapshot()['cycle_sample_count'],0)
        time.sleep(.4)
        with self.log.open('a') as f:f.write('Step number 3 out of a maximum of 40\n')
        m.poll()
        self.assertEqual(m.snapshot()['cycle_sample_count'],1)
        self.assertGreater(m.snapshot()['cycle_average_s'],0)

    def test_completed_attached_log_shows_no_future_prediction(self):
        self.log.write_text(''' Leave Link 103 at Sun Oct 18 04:01:00 2026, cpu: 1.0 elap: 1.0
 Step number 1 out of a maximum of 30
 Leave Link 103 at Sun Oct 18 04:01:20 2026, cpu: 1.0 elap: 1.0
 Step number 2 out of a maximum of 30
 Normal termination of Gaussian 16 at Sun Oct 18 04:02:00 2026.
''')
        m=Monitor();m.attach(self.log)
        s=m.poll()
        self.assertEqual(s['state'],'COMPLETED')
        self.assertIsNone(s['eta_records'][0]['seconds'])
