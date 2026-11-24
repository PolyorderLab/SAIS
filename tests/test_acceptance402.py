"""Acceptance checks: genuine motifs and contradictory geometry/spectra."""
import sys,unittest
from pathlib import Path
import numpy as np
sys.path[:0]=[str(Path(__file__).resolve().parents[1]/'structure061')]
from acceptance402 import apply,PROTOTYPES,periodic_euler,network_check,prototype_fit

def spheres(points,shape=(36,36,36),A=None,sigma=.10,weights=None):
    A=np.eye(3) if A is None else np.asarray(A)
    coords=np.stack(np.meshgrid(*[np.arange(n)/n for n in shape],indexing='ij'),axis=-1)
    v=np.zeros(shape)
    for p,w in zip(points,weights if weights is not None else np.ones(len(points))):
        d=(coords-p+.5)%1-.5;v+=w*np.exp(-np.sum((d@A.T)**2,axis=-1)/(2*sigma**2))
    v=.01+.95*v/v.max()
    return [v,1-v],A

class AcceptanceTests(unittest.TestCase):
    def accepted(self,name,f,A):return apply(dict(candidate=name,status='scattering_supported'),f,A)
    def test_simple_cubic_and_origin(self):
        for shift in ([0,0,0],[.25,.5,.75]):
            f,A=spheres(np.asarray(PROTOTYPES['SC'])+shift)
            self.assertEqual(self.accepted('SC',f,A)['candidate'],'SC')
    def test_body_and_face_centered(self):
        for name in ['BCC','FCC']:
            f,A=spheres(PROTOTYPES[name],sigma=.065)
            self.assertEqual(self.accepted(name,f,A)['candidate'],name)
    def test_sphere_forming_majority_is_not_automatically_rejected(self):
        f,A=spheres(PROTOTYPES['BCC'],sigma=.22)
        # A constant offset changes which component is the majority, not its
        # relative-level morphology or nonzero Fourier modes.
        f=[f[0]+.4, f[1]-.4]
        self.assertEqual(self.accepted('BCC',f,A)['candidate'],'BCC')
    def test_extra_satellites_cannot_pass_fcc(self):
        p=np.array(list(__import__('itertools').product([0,.5],repeat=3)))
        w=[1 if list(x) in PROTOTYPES['FCC'].tolist() else .4 for x in p]
        f,A=spheres(p,sigma=.06,weights=w)
        self.assertIsNone(self.accepted('FCC',f,A)['candidate'])
    def test_unequal_fcc_sites_cannot_pass_centering(self):
        f,A=spheres(PROTOTYPES['FCC'],sigma=.065,weights=[1,1,1,.98])
        self.assertIsNone(self.accepted('FCC',f,A)['candidate'])
    def test_anisotropic_ball_not_standard_sc(self):
        f,A=spheres(PROTOTYPES['SC']);z=np.arange(36)/36
        f=[f[0]+.03*np.cos(2*np.pi*z)[:,None,None],f[1]-.03*np.cos(2*np.pi*z)[:,None,None]]
        self.assertIsNone(self.accepted('SC',f,A)['candidate'])
    def test_non_cubic_motifs_and_metric(self):
        for name,A in [('HCP',np.diag([1,np.sqrt(3),2.4])),('iHPa',np.diag([1,np.sqrt(3),1.1])),('PtS',np.diag([1,1,3.0]))]:
            self.assertLess(prototype_fit(PROTOTYPES[name]+.13,A,name),1e-10)
        self.assertIsNone(prototype_fit(PROTOTYPES['HCP'],np.eye(3),'HCP'))
    def test_lam_diagnostic_failure_is_not_bypassed(self):
        r=apply(dict(candidate='LAM',status='one_dimensional_structure_supported',lam_diagnostic=dict(family_supported=False)),[np.arange(16),-np.arange(16)],np.eye(1))
        self.assertIsNone(r['candidate'])
    def test_review_is_not_promoted(self):
        f,A=spheres(PROTOTYPES['BCC'])
        self.assertEqual(apply(dict(candidate='BCC',status='review_reflection_or_threshold_conflict'),f,A)['status'],'review_reflection_or_threshold_conflict')
    def test_periodic_euler_is_translation_invariant(self):
        f,A=spheres(PROTOTYPES['SC']);mask=f[0]>.5
        self.assertEqual(periodic_euler(mask),1)
        self.assertEqual(periodic_euler(np.roll(mask,13,axis=1)),1)
        slab=np.zeros((20,20,20),bool);slab[:5]=True
        self.assertEqual(periodic_euler(slab),0)
    def test_spheres_are_not_gyroid_or_diamond_networks(self):
        f,A=spheres(PROTOTYPES['FCC'])
        for name in ('SG','DG','DD'):self.assertFalse(network_check(f,A,name)['passed'])
    def test_single_gyroid_necessary_signature(self):
        x,y,z=np.meshgrid(*[np.arange(40)*2*np.pi/40]*3,indexing='ij')
        g=np.sin(x)*np.cos(y)+np.sin(y)*np.cos(z)+np.sin(z)*np.cos(x)
        v=(g-g.min())/np.ptp(g)
        self.assertTrue(network_check([v,1-v],np.eye(3),'SG')['passed'])

if __name__=='__main__':unittest.main()
