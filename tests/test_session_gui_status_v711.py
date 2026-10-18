"""Regression for window crash with saved-but-unattached job history."""
from dataclasses import dataclass
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from gaumon_ext.session import SessionStore,describe_observed_job

class SessionStatusTests(unittest.TestCase):
    def test_saved_record_has_no_loaded_snapshot(self):
        self.assertEqual(describe_observed_job(None),'SAVED · RESUME TO ATTACH')

    def test_uninitialized_and_unreadable_log_reader(self):
        x=SimpleNamespace(snapshot=None,last_poll_error='',last_growth_mono=None)
        self.assertEqual(describe_observed_job(x),'ATTACHING · WAITING FOR LOG')
        x.last_poll_error='Permission denied'
        self.assertEqual(describe_observed_job(x),'UNREADABLE')

    def test_unknown_vs_completed_vs_failed(self):
        x=SimpleNamespace(snapshot=SimpleNamespace(termination=''),last_poll_error='',last_growth_mono=None)
        self.assertEqual(describe_observed_job(x),'LOG ATTACHED · PROCESS UNVERIFIED')
        x.snapshot.termination='FAILED'
        self.assertEqual(describe_observed_job(x),'FAILED')
        x.snapshot.termination='COMPLETED'
        self.assertEqual(describe_observed_job(x),'COMPLETED')

    def test_recent_log_growth_is_not_process_proof(self):
        x=SimpleNamespace(snapshot=SimpleNamespace(termination=''),last_poll_error='',last_growth_mono=97.)
        self.assertEqual(describe_observed_job(x,now=105.),'LOG GROWING · PROCESS UNVERIFIED')
        self.assertEqual(describe_observed_job(x,now=128.),'LOG ATTACHED · PROCESS UNVERIFIED')

    def test_saved_job_restores_path_not_active_log_reader(self):
        with tempfile.TemporaryDirectory() as tmp:
            store_file=Path(tmp)/'sessions.json'
            log=Path(tmp)/'a.log';log.write_text('Entering Gaussian System\n')
            store=SessionStore(store_file);store.add(log)
            loaded=SessionStore(store_file)
            self.assertEqual(len(loaded.items),1)
            self.assertEqual(describe_observed_job(None),'SAVED · RESUME TO ATTACH')
            self.assertIsNone(loaded.items[0]['pid'])

if __name__=='__main__':unittest.main()
