"""Scientific parsing regression tests with synthetic data, no licensed Gaussian."""
import math
import unittest
from pathlib import Path
from gaumon.engine import Parser, Monitor
from gaumon.deep_scf import DeepSCFParser, MAX_ROWS_PER_BLOCK
from gaumon.inputs import ADVANCED_CASES, ALL_CASES, load_case, validate_input

ROOT=Path(__file__).resolve().parents[1]
STRESS=ROOT/'examples'/'stress_logs'


class TestDeepParsing(unittest.TestCase):
    def parse(self,text):
        parser=DeepSCFParser()
        for line in text.splitlines():parser.feed_line(line)
        return parser

    def test_user_d_notation_values(self):
        text=''' Cycle   7  Pass 1  IDiag  1:
 RMSU= 9.47D-04 CP: 9.64D-01 1.03D-02 7.27D-01 3.97D-01 3.16D-01
                    CP: 8.35D-01
 E= -4599.34852544226 Delta-E= -0.562219918020 Rises=F Damp=F
 DIIS: error= 1.03D-02 at cycle 7 NSaved= 7.
 NSaved= 7 IEnMin= 7 EnMin= -4599.34852544226 IErMin= 7 ErrMin= 1.03D-02
 IDIUse=3 WtCom= 8.97D-01 WtEn= 1.03D-01
 Coeff-Com: -0.165D+00 0.225D+00 0.274D-01-0.804D-01-0.321D-01 0.257D+00
 Coeff-Com: 0.767D+00
 Coeff-En: 0.000D+00 0.000D+00 0.100D+01
 Coeff: -0.148D+00 0.202D+00 0.246D-01-0.721D-01
 Gap= 0.330 Goal= None Shift= 0.000
 Gap= 0.286 Goal= None Shift= 0.000
 RMSDP=5.13D-04 MaxDP=5.66D-02 DE=-5.62D-01 OVMax= 5.82D-02
 Cycle   8  Pass 1 IDiag 1:
 E= -4599.43187844636 Delta-E= -0.083353004100 Rises=F Damp=F
 DIIS: error= 4.15D-03 at cycle 8 NSaved= 8.
 RMSDP=5.03D-04 MaxDP=4.45D-02 DE=-8.34D-02 OVMax= 7.75D-02
'''
        p=self.parse(text)
        b=p.get_block()
        self.assertEqual(len(b['rows']),2)
        r=b['rows'][0]
        self.assertAlmostEqual(r['energy'],-4599.34852544226,10)
        self.assertAlmostEqual(r['delta_e'],-0.562219918020,10)
        self.assertAlmostEqual(r['diis_error'],0.0103)
        self.assertAlmostEqual(r['rmsu'],0.000947)
        self.assertAlmostEqual(r['rmsdp'],0.000513)
        self.assertAlmostEqual(r['maxdp'],0.0566)
        self.assertAlmostEqual(r['ovmax'],0.0582)
        self.assertEqual(r['gaps'],[0.330,0.286])
        self.assertEqual(len(r['cp']),6)
        self.assertEqual(len(r['coeff_com']),7)
        self.assertAlmostEqual(r['coeff_com'][3],-0.0804)
        self.assertEqual(r['nsaved'],7.0)
        self.assertEqual(b['rows'][1]['cycle'],8)

    def test_incremental_one_character_chunks_equivalent(self):
        text=(STRESS/'chunked.log').read_text()
        a=Parser();a.feed(text)
        b=Parser()
        for n in range(0,len(text),3):b.feed(text[n:n+3])
        self.assertEqual(a.snapshot()['deep_scf'], b.snapshot()['deep_scf'])
        self.assertEqual(a.deep.get_block(1), b.deep.get_block(1))

    def test_reset_and_history(self):
        text=(STRESS/'many_steps.log').read_text()
        parser=Parser();parser.feed(text)
        m=parser.snapshot()['deep_scf']
        self.assertEqual(len(m['blocks']),8)
        self.assertTrue(all(b['status']=='complete' for b in m['blocks']))
        self.assertEqual([b['optimization_step'] for b in m['blocks']],list(range(1,9)))
        self.assertEqual(m['active_block_id'],8)
        self.assertEqual(m['active_cycle'],None)
        self.assertEqual(parser.deep.get_block(1)['rows'][0]['cycle'],1)
        self.assertEqual(parser.deep.get_block(8)['rows'][0]['cycle'],1)
        self.assertEqual(len(parser.deep.get_block(8)['rows']),17)

    def test_oscillating_convergence_stress(self):
        parser=Parser();parser.feed((STRESS/'oscillation.log').read_text())
        m=parser.snapshot()['deep_scf']
        self.assertEqual(len(m['blocks']),1)
        self.assertEqual(m['active_cycle'],95)
        block=parser.deep.get_block()
        self.assertEqual(block['status'],'running')
        self.assertTrue(any(r.get('rises')=='T' for r in block['rows']))
        self.assertTrue(all(math.isfinite(r['energy']) for r in block['rows']))
        self.assertIsNone(block['final_energy'])

    def test_abrupt_eof_partial_line_no_false_completion(self):
        parser=Parser();parser.feed((STRESS/'incomplete.log').read_text())
        block=parser.deep.get_block()
        self.assertEqual(block['status'],'running')
        self.assertEqual(len(block['rows']),12)
        self.assertIsNone(block['final_energy'])

    def test_partial_next_step(self):
        parser=Parser();parser.feed((STRESS/'chunked.log').read_text())
        m=parser.snapshot()['deep_scf']
        self.assertEqual(len(m['blocks']),2)
        self.assertEqual(m['active_cycle'],1)
        self.assertEqual(parser.deep.get_block(1)['optimization_step'],1)
        self.assertEqual(parser.deep.get_block(2)['status'],'running')

    def test_repeated_header_same_pass_is_not_duplicate(self):
        p=self.parse(' Cycle 1 Pass 1 IDiag 1:\n E= -2\n Cycle 1 Pass 1 IDiag 1:\n DIIS: error=1.00D-03 at cycle 1\n')
        self.assertEqual(len(p.get_block()['rows']),1)
        self.assertAlmostEqual(p.get_block()['rows'][0]['diis_error'],.001)

    def test_max_rows_bounded(self):
        parser=DeepSCFParser()
        for i in range(1,MAX_ROWS_PER_BLOCK+80):parser.feed_line(f' Cycle {i} Pass 1 IDiag 1:')
        block=parser.get_block()
        self.assertEqual(len(block['rows']),MAX_ROWS_PER_BLOCK)
        self.assertTrue(block['truncated'])
        self.assertEqual(block['rows'][-1]['cycle'],MAX_ROWS_PER_BLOCK+79)

    def test_snapshots_do_not_leak_mutable_data(self):
        parser=self.parse(' Cycle 1 Pass 1 IDiag 1:\n E= -5.0\n')
        block=parser.get_block();block['rows'][0]['energy']=123.0
        self.assertEqual(parser.get_block()['rows'][0]['energy'],-5.)

    def test_advanced_gaussian_inputs_valid_but_not_benchmarked(self):
        self.assertEqual(len(ALL_CASES),7)
        self.assertEqual(len(ADVANCED_CASES),3)
        for case in ADVANCED_CASES:
            with self.subTest(case=case.key):
                self.assertFalse(validate_input(load_case(case.key)[1]))
                self.assertNotIn('geom=connectivity', load_case(case.key)[1].lower())


if __name__=='__main__':unittest.main()
