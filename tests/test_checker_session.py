import tempfile,unittest
from pathlib import Path
from gaumon_ext.checker import check_input
from gaumon_ext.session import SessionStore
class CheckTests(unittest.TestCase):
 def test_valid_gaussian_connectivity(self):
  text='''%mem=1GB\n#p hf/3-21g opt geom=connectivity\n\nWater\n\n0 1\nO 0 0 0\nH 0 0 1\n\n1 2 1.0\n2\n\n'''
  problems=check_input(text)
  self.assertEqual([p for p in problems if p.severity=='error'],[])
 def test_missing_separator(self):
  text='#p hf/3-21g opt\nWater\n0 1\nH 0 0 0\n'
  issues=check_input(text)
  self.assertTrue(any(p.severity=='error' for p in issues))
 def test_trailing_spaces_visible(self):
  issues=check_input('#p hf/sto-3g sp   \n\nTest\n\n0 1\nH 0 0 0\n\n')
  self.assertTrue(any('Trailing' in p.message for p in issues))
 def test_link1_sections(self):
  text='# hf/3-21g sp\n\nFirst\n\n0 1\nH 0 0 0\n\n--Link1--\n# hf/3-21g sp\n\nSecond\n\n0 1\nH 0 0 0\n\n'
  self.assertFalse(any(p.severity=='error' for p in check_input(text)))
 def test_session_persists_only_explicit_paths(self):
  with tempfile.TemporaryDirectory() as d:
   f=Path(d)/'log.log';f.write_text('a')
   db=Path(d)/'history.json'
   a=SessionStore(db);a.add(f);a.add(f)
   self.assertEqual(len(a.items),1)
   b=SessionStore(db);self.assertEqual(len(b.items),1)
   self.assertIsNone(b.items[0]['pid'])
   b.remove(b.items[0]['id']);self.assertEqual(SessionStore(db).items,[])
