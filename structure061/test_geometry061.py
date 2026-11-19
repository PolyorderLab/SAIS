import unittest
import numpy as np
from geometry061 import lamellar,planar

class GeometryTests(unittest.TestCase):
 def test_fractional_grid_period(self):
  a=.3+.2*np.cos(4*np.pi*np.arange(135)/135)
  self.assertIsNotNone(lamellar([a,1-a]))
 def test_unequal_layers_rejected(self):
  x=np.arange(135)/135;a=.3+.15*np.cos(4*np.pi*x)+.08*np.cos(2*np.pi*x)
  self.assertIsNone(lamellar([a,1-a]))
 def test_diagonal_and_two_direction_fields(self):
  u,v=np.meshgrid(np.arange(64)/64,np.arange(64)/64,indexing='ij')
  a=.3+.2*np.cos(2*np.pi*(u+v));self.assertIsNotNone(lamellar([a,1-a]))
  a=.3+.1*np.cos(2*np.pi*u)+.1*np.cos(2*np.pi*v);self.assertIsNone(lamellar([a,1-a]))
 def test_ring_rejected(self):
  u,v=np.meshgrid(np.arange(64)/64,np.arange(64)/64,indexing='ij');A=np.array([[1,-.5],[0,np.sqrt(3)/2]])
  rr=np.linalg.norm(np.stack([u-.5,v-.5],axis=-1)@A.T,axis=-1)
  a=.2+.6*np.exp(-((rr-.2)/.04)**2);self.assertFalse(planar([a,1-a],A))
 def test_real_wrapped_HEX_and_ring(self):
  from pathlib import Path
  from run import load
  p=Path(__file__).resolve().parent.parent/'validation057/inputs'
  for i in [6702,6711,7428]:
   fs,A,_=load(p/f'run{i}');self.assertTrue(planar(fs,A))
  fs,A,_=load(p/'run7272');self.assertFalse(planar(fs,A))

if __name__=='__main__':unittest.main()
