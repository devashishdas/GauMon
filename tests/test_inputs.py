import tempfile
import unittest
from pathlib import Path
from gaumon.inputs import CASES,load_case,write_case,validate_input
from gaumon.diagnostics import diagnose_log

class TestGaussianInputs(unittest.TestCase):
    def test_all_four_templates_are_well_formed(self):
        self.assertEqual(len(CASES),4)
        for c in CASES:
            with self.subTest(c=c.key):
                case,text=load_case(c.key)
                self.assertEqual(validate_input(text),[])
                self.assertTrue(text.count('\n\n')>=2)
                self.assertIn('GauMon test',text)
                self.assertNotIn('geom=connectivity',text)

    def test_missing_title_from_user_error_is_detected(self):
        malformed='#p opt freq b3lyp/aug-cc-pvtz geom=connectivity\n\n\n0 1\nO 0 0 0\nH 0 0 1\n\n'
        issues=validate_input(malformed)
        self.assertTrue(any('Title Card Required' in i for i in issues),issues)

    def test_missing_charge_or_coordinates(self):
        self.assertTrue(any('charge/multiplicity' in i for i in validate_input('#p hf/3-21g\n\nExample\n\nH 0 0 0\n')))
        self.assertTrue(any('coordinates' in i for i in validate_input('#p hf/3-21g\n\nExample\n\n0 1\n')))

    def test_missing_connectivity_block(self):
        text='#p hf/3-21g geom=connectivity\n\nExample\n\n0 1\nH 0 0 0\nH 0 0 1\n\n'
        self.assertTrue(any('connectivity' in i for i in validate_input(text)))

    def test_case_export_no_silent_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=write_case('sp',Path(tmp))
            self.assertTrue(p.is_file())
            with self.assertRaises(FileExistsError):write_case('sp',Path(tmp))
            p2=write_case('sp',Path(tmp),overwrite=True)
            self.assertEqual(p,p2)

    def test_error_is_prior_to_segfault(self):
        data='(Enter /g16/l101.exe)\nTitle Card Required\nEnd of file in ZSymb.\nError termination via Lnk1e in l101.exe\nError: segmentation violation\n'
        messages=diagnose_log(data)
        self.assertIn('ZSymb',messages[0].title)
        self.assertEqual(messages[1].severity,'INFO')
        self.assertTrue(any(m.title=='Segmentation fault reported' for m in messages))

    def test_missing_final_blank_line_is_detected(self):
        content=load_case('sp')[1].rstrip('\n')+'\n'
        self.assertTrue(any('ZSymb risk' in issue for issue in validate_input(content)))

    def test_title_card_required_alone_is_not_a_fatal_error(self):
        diagnosis=diagnose_log('Title Card Required\nNormal termination of Gaussian 16\n')
        self.assertTrue(diagnosis)
        self.assertEqual(diagnosis[0].severity,'INFO')

    def test_valid_important_case_routes(self):
        self.assertIn('opt freq',load_case('optfreq')[1].lower())
        self.assertIn('opt=(ts,calcfc,noeigentest',load_case('ts')[1].lower())
        self.assertIn('3-21g freq', load_case('freq')[1].lower())
