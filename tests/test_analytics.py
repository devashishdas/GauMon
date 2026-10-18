import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from gaumon_ext.analytics import IncrementalLog,inspect_text

HERE=Path(__file__).resolve().parents[1]
class AnalyticsTests(unittest.TestCase):
 def test_large_link_history_is_bounded_but_totals_complete(self):
  header='Leave Link  502 at Sun Oct 18 11:00:18 2026, MaxMem= 123 cpu: 2.0 elap: 1.0\n'
  snap=inspect_text(header*2200)
  self.assertEqual(len(snap.events),2000)
  self.assertEqual(snap.event_count,2200)
  self.assertEqual(snap.totals()[502],[2200.0,4400.0,2200])
 def test_link_times_and_errors(self):
  s=inspect_text((HERE/'examples/mock_timing.log').read_text())
  self.assertEqual(len(s.events),4)
  self.assertEqual(s.current_link,103)
  self.assertAlmostEqual(s.totals()[502][0],17)
  self.assertAlmostEqual(s.totals()[502][1],55)
  self.assertEqual(s.totals()[502][2],2)
  self.assertEqual(s.termination,'FAILED')
  self.assertTrue(any('FormBX' in x.evidence for x in s.problems))
  self.assertFalse(s.eta()['available'])
 def test_incremental_partial_chunks(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'test.log';p.write_bytes(b'(Enter /x/l502.exe)\nLeave Link  502 at Sun Oct 18 11:00:18 2026, MaxMem= 123 cpu:')
   engine=IncrementalLog(p);engine.poll(max_bytes=31)
   self.assertEqual(len(engine.snapshot.events),0)
   while engine.offset<p.stat().st_size:engine.poll(max_bytes=31)
   with p.open('ab') as out:out.write(b' 21.0 elap: 7.0\n')
   engine.poll(max_bytes=31)
   self.assertEqual(engine.snapshot.totals()[502],[7.0,21.0,1])
 def test_truncate_resets_records(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'x.log';p.write_text((HERE/'examples/mock_timing.log').read_text())
   e=IncrementalLog(p)
   while e.offset<p.stat().st_size:e.poll()
   self.assertEqual(len(e.snapshot.events),4)
   p.write_text('tiny\n');e.poll()
   self.assertEqual(len(e.snapshot.events),0)
 def test_eta_abstains_on_backfilled_log(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'x.log';p.write_text('Step number 1 out of a maximum of 50\nStep number 2 out of a maximum of 50\n')
   e=IncrementalLog(p);e.poll()
   self.assertEqual(list(e.snapshot.observed_step_times),[])
   self.assertFalse(e.snapshot.eta()['available'])
 def test_live_eta_is_provisional(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'x.log';p.write_text('header\n')
   e=IncrementalLog(p);e.poll()
   values=iter([10.0,10.0,20.0,20.0,32.0,32.0])
   with patch('gaumon_ext.analytics.time.monotonic',side_effect=lambda:next(values)):
    for i in (1,2,3):
     with p.open('a') as out:out.write(f'Step number {i} out of a maximum of 60\n')
     e.poll()
   self.assertTrue(e.snapshot.eta()['available'])
   self.assertEqual(e.snapshot.eta()['sample_count'],2)
   self.assertAlmostEqual(e.snapshot.eta()['pace_s'],11.)
