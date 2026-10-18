"""No-Qt regression coverage for discovery, route help, deque ETA, and persistence."""
import tempfile
import unittest
from collections import deque
from pathlib import Path
from gaumon_ext.analytics import adjacent_step_durations, Snapshot
from gaumon_ext.input_context import matching_input_for_log,explain_input,keyword_info
from gaumon_ext.session import SessionStore

class InputContextTests(unittest.TestCase):
    def test_deque_iteration_regression_repeated_updates(self):
        observations=deque([(1,1.),(2,4.),(3,8.),(8,10.),(9,19.)],maxlen=500)
        for repeat in range(20):
            self.assertEqual(adjacent_step_durations(observations),[(2,3.),(3,4.),(9,9.)])
        s=Snapshot();s.observed_step_times=observations
        self.assertTrue(s.eta()['available'])
        self.assertEqual(s.eta()['sample_count'],3)

    def test_explain_gaussian_route_and_link0(self):
        inp='%nprocshared=4\n%mem=8GB\n#p ub3lyp/def2svp sp scf=(xqc,maxcycle=192) nosymm\n\nFe2S2\n\n-2 1\nFe 0 0 0\n\n'
        records=explain_input(inp)
        names=[r[1].lower() for r in records]
        self.assertIn('%nprocshared',names)
        self.assertIn('%mem',names)
        self.assertIn('#p',names)
        for v in ['ub3lyp','def2svp','sp','scf','xqc','maxcycle','nosymm']:self.assertIn(v,names)
        self.assertFalse(any(x[1]=='Fe2S2' for x in records))
        self.assertIsNone(keyword_info('arbitraryexperimentalmethod'))

    def test_separate_link1_sections(self):
        inp='#p sp hf/3-21g\n\ntest\n\n0 1\nH 0 0 0\n\n--Link1--\n#p opt freq\n\ntest\n\n0 1\nH 0 0 0\n\n'
        names=[r[1].lower() for r in explain_input(inp)]
        self.assertEqual(names.count('sp'),1)
        self.assertEqual(names.count('freq'),1)
        self.assertEqual(names.count('opt'),1)

    def test_only_matching_original_not_arbitrary_gjf(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);(d/'remote').mkdir()
            log=d/'remote'/'fe2s2.log';log.write_text('no source path printed\n')
            (d/'fe2s2.gjf').write_text('wrong file\n')
            self.assertIsNone(matching_input_for_log(log,cwd=d))
            correct=log.with_suffix('.gjf');correct.write_text('correct\n')
            self.assertEqual(matching_input_for_log(log,cwd=d),correct)
            (d/'remote'/'x.log.txt').write_text('test\n')
            x=d/'remote'/'x.com';x.write_text('correct\n')
            self.assertEqual(matching_input_for_log(d/'remote'/'x.log.txt'),x)

    def test_explicit_session_input_provenance_survives_restart(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);log=d/'x.log';log.write_text('')
            gjf=d/'other.gjf';gjf.write_text('#p hf/sto-3g sp\n\ntitle\n\n0 1\nH 0 0 0\n\n')
            st=SessionStore(d/'state.json');record=st.add(log,'started-in-GauMon',gjf)
            self.assertEqual(record['input_path'],str(gjf.resolve()))
            again=SessionStore(d/'state.json')
            self.assertEqual(again.items[0]['input_path'],str(gjf.resolve()))
            self.assertEqual(again.add(log)['id'],record['id'])
            self.assertEqual(again.items[0]['input_path'],str(gjf.resolve()))

if __name__=='__main__':unittest.main()
