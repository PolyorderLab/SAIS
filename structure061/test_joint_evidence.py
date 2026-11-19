import unittest,tempfile,json
from pathlib import Path
import numpy as np
import yaml
from joint_evidence import apply

class EvidenceGateTests(unittest.TestCase):
 def fixture(self,amp=.003,F=7.64000013,status='JuPolyorder.Successful()',residual=1e-7):
  import uuid
  p=Path(__file__).resolve().parent/'test_artifacts'/uuid.uuid4().hex;p.mkdir(parents=True)
  cfg={'system':{'components':[{'type':'BCP','volume_fraction':1.,'length':1.,'blocks':[{'segment':'A','length':.2},{'segment':'B','length':.8}]}],'chain_density':1.,'χN_map':[['A','B',54.]]},'scft':{'tolmode':'Residual','relative':False,'tol':1e-4},'cellopt':{'gtol':1e-5}}
  summary={'F':F,'residual':residual,'stress_norm':1e-6,'convergence':status}
  (p/'config.yml').write_text(yaml.safe_dump(cfg,allow_unicode=True),encoding="utf8");(p/'summary.yml').write_text(yaml.safe_dump(summary))
  a=.2+amp*np.cos(2*np.pi*np.arange(80)/40);return p,[a,1-a],cfg
 def test_weak_wave_near_DIS_not_forced_to_either_phase(self):
  p,fs,_=self.fixture();before=[a.copy() for a in fs];r=apply({'candidate':'LAM','status':'one_dimensional_structure_supported'},fs,p)
  self.assertEqual(r['status'],'review_DIS_vs_weak_order');self.assertIsNone(r['candidate']);self.assertEqual(r['morphology_candidate'],'LAM')
  for a,b in zip(fs,before):np.testing.assert_array_equal(a,b)
 def test_uniform_energy_agreement(self):
  p,fs,_=self.fixture(amp=1e-9,F=7.64);r=apply({'candidate':None,'status':'DIS_candidate'},fs,p);self.assertEqual(r['candidate'],'DIS')
 def test_failed_convergence_blocks_strong_order(self):
  p,fs,_=self.fixture(amp=.15,F=7.,status='JuPolyorder.Failed()');r=apply({'candidate':'LAM','status':'one_dimensional_structure_supported'},fs,p);self.assertEqual(r['status'],'review_convergence_unverified')
 def test_reported_success_with_large_residual_is_not_enough(self):
  p,fs,_=self.fixture(amp=.15,F=7.,residual=.1);self.assertEqual(apply({'candidate':'LAM'},fs,p)['status'],'review_convergence_unverified')
 def test_lower_energy_does_not_override_scattering_conflict(self):
  p,fs,_=self.fixture(amp=.15,F=7.);r=apply({'candidate':'HEX','status':'review_reflection_or_threshold_conflict'},fs,p);self.assertIsNone(r['candidate'])
 def test_strong_order_below_DIS_is_proposal_not_stability_proof(self):
  p,fs,_=self.fixture(amp=.15,F=7.);r=apply({'candidate':'LAM','status':'one_dimensional_structure_supported'},fs,p);self.assertEqual(r['candidate'],'LAM');self.assertFalse(r['joint_evidence']['thermodynamic_global_stability_proven']);self.assertFalse(r['database_written'])
 def test_weak_order_below_DIS_not_misclassified_DIS(self):
  p,fs,_=self.fixture(amp=.003,F=7.63);r=apply({'candidate':'LAM'},fs,p);self.assertEqual(r['status'],'review_weak_order_below_DIS');self.assertNotEqual(r['candidate'],'DIS')
 def test_missing_config_preserves_morphology_only(self):
  p,fs,_=self.fixture();(p/'config.yml').unlink();r=apply({'candidate':'LAM'},fs,p);self.assertEqual(r['morphology_candidate'],'LAM');self.assertIsNone(r['candidate'])
 def test_config_density_mismatch_rejects_baseline(self):
  p,fs,cfg=self.fixture();cfg['system']['components'][0]['blocks']=[{'segment':'A','length':.3},{'segment':'B','length':.7}];(p/'config.yml').write_text(yaml.safe_dump(cfg));r=apply({'candidate':'LAM'},fs,p);self.assertFalse(r['joint_evidence']['energy']['available'])
 def test_above_DIS_does_not_prove_DIS(self):
  p,fs,_=self.fixture(amp=.15,F=7.7);r=apply({'candidate':'LAM'},fs,p);self.assertEqual(r['status'],'review_above_DIS_energy');self.assertIsNone(r['candidate'])
if __name__=='__main__':unittest.main()
