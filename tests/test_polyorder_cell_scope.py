import copy,json,math,sys,tempfile,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'investigation068'))
from base061 import convergence_evidence,apply
class CellScopeTests(unittest.TestCase):
 def config(self):
  return {'scft':{'tol':1e-7,'relative':False,'tolmode':'Residual'},'cellopt':{'tol_stress':1e-5},'system':{'components':[{'type':'BCP','volume_fraction':1.,'length':1.,'blocks':[{'segment':'A','length':.4},{'segment':'B','length':.6}]}],'chain_density':1.,'χN_map':[['A','B',30.]]}}
 def summary(self):return {'convergence':'Polyorder.Successful()','residual':1e-8,'stress_norm':.3,'F':5.}
 def test_fixed_large_stress_is_not_stationarity_test(self):
  d=convergence_evidence(self.summary(),self.config(),{'kind':'fixed_cell'})
  self.assertTrue(d['passed']);self.assertFalse(d['stress_required']);self.assertEqual(d['stress_norm'],.3)
 def test_optimized_large_stress_rejected(self):
  self.assertFalse(convergence_evidence(self.summary(),self.config(),{'kind':'cell_optimized'})['passed'])
 def test_optimized_passes_only_both(self):
  s=self.summary();s['stress_norm']=1e-6
  self.assertTrue(convergence_evidence(s,self.config(),{'kind':'cell_optimized'})['passed'])
  s['residual']=1e-6;self.assertFalse(convergence_evidence(s,self.config(),{'kind':'cell_optimized'})['passed'])
 def test_no_scope_remains_strict(self):
  self.assertFalse(convergence_evidence(self.summary(),self.config())['passed'])
 def test_legacy_success_and_gtol(self):
  c=self.config();c['cellopt']={'gtol':1e-5};s=self.summary();s.update(convergence='JuPolyorder.Successful()',stress_norm=1e-6)
  self.assertTrue(convergence_evidence(s,c)['passed'])
 def test_unknown_or_invalid_scope_rejected(self):
  for context in [{'kind':'unknown'},{},[],{'kind':None}]:
   self.assertFalse(convergence_evidence(self.summary(),self.config(),context)['passed'])
 def test_success_status_not_loosened(self):
  for status in ['Polyorder.Acceptable()','JuPolyorder.Acceptable()','Successful','Polyorder.NeedMoreIteration()',None]:
   s=self.summary();s['convergence']=status
   self.assertFalse(convergence_evidence(s,self.config(),{'kind':'fixed_cell'})['passed'])
 def test_nonfinite_negative_or_excessive_residual(self):
  for value in [-1.,float('nan'),float('inf'),1e-5]:
   s=self.summary();s['residual']=value;self.assertFalse(convergence_evidence(s,self.config(),{'kind':'fixed_cell'})['passed'])
 def test_relative_stop_not_accepted(self):
  c=self.config();c['scft']['relative']=True
  self.assertFalse(convergence_evidence(self.summary(),c,{'kind':'fixed_cell'})['passed'])
 def test_conflicting_aliases_rejected(self):
  c=self.config();c['cellopt']['gtol']=1e-3;s=self.summary();s['stress_norm']=1e-6
  self.assertFalse(convergence_evidence(s,c,{'kind':'cell_optimized'})['passed'])
 def test_missing_stress_only_ok_fixed(self):
  s=self.summary();s.pop('stress_norm')
  self.assertTrue(convergence_evidence(s,self.config(),{'kind':'fixed_cell'})['passed'])
  self.assertFalse(convergence_evidence(s,self.config(),{'kind':'cell_optimized'})['passed'])
 def test_json_input_and_other_gates_unchanged(self):
  v=.4+.15*np.cos(np.arange(64)*2*np.pi/64)
  with tempfile.TemporaryDirectory() as temp:
   d=Path(temp);c=self.config();s=self.summary()
   (d/'config.json').write_text(json.dumps(c));(d/'summary.json').write_text(json.dumps(s));(d/'calculation.json').write_text(json.dumps({'kind':'fixed_cell'}))
   result=apply({'candidate':'LAM','status':'one_dimensional_structure_supported'},[v,1-v],d)
   self.assertEqual(result['candidate'],'LAM');self.assertEqual(result['status'],'joint_ordered_supported')
   self.assertIn(str(d/'calculation.json'),result['joint_evidence']['metadata_hashes'])
   s['F']=7.;(d/'summary.json').write_text(json.dumps(s))
   result=apply({'candidate':'LAM','status':'one_dimensional_structure_supported'},[v,1-v],d)
   self.assertEqual(result['status'],'review_above_DIS_energy');self.assertIsNone(result['candidate'])
   s['F']=5.;(d/'summary.json').write_text(json.dumps(s));(d/'calculation.json').write_text('[]')
   self.assertIsNone(apply({'candidate':'LAM','status':'one_dimensional_structure_supported'},[v,1-v],d)['candidate'])
if __name__=='__main__':unittest.main()
