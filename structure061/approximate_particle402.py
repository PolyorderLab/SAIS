"""Restricted approximate axis-supercell route for an existing SC scattering hypothesis.
No file names, seed labels or sample identifiers participate.
"""
import itertools
import numpy as np
from scipy.optimize import linear_sum_assignment

POLICY=dict(axis_spectrum_proposal_threshold=1e-8,maximum_discarded_power_fraction=1e-8,
 maximum_normalized_point_error=1e-3,minimum_primitive_grid_points=12,
 maximum_domains=64,relative_levels=[.25,.35,.5,.65,.75],
 center_tolerance_primitive_fraction=.02,maximum_size_cv=.15,maximum_axis_ratio=4.)

def evaluate(fields,A):
 from supercell402 import reduce_supercell
 from particle402 import domains,describe
 a=[np.asarray(x) for x in fields];A=np.asarray(A)
 ev=dict(policy=POLICY,accepted=False,scope='approximate axis-repeat SC particle hypothesis; not exact symmetry or global stability proof')
 if len(a)!=2 or A.shape!=(3,3) or any(x.ndim!=3 or not np.isfinite(x).all() for x in a):return None,ev
 small,cell,validation=reduce_supercell(a,A,power_threshold=POLICY['axis_spectrum_proposal_threshold'],relative_tolerance=POLICY['maximum_normalized_point_error'])
 ev['translation_validation']=validation
 if not validation.get('reduced'):ev['reason']='no_approximate_repeat_proposal';return None,ev
 repeats=np.array(validation['repeats'],dtype=int);n=int(np.prod(repeats));shape=np.array(a[0].shape)
 if n<=1 or n>POLICY['maximum_domains'] or min(shape/repeats)<POLICY['minimum_primitive_grid_points']:
  ev['reason']='repeat_count_or_primitive_sampling_unsupported';return None,ev
 lengths=np.linalg.norm(cell,axis=0)
 metric=max(np.ptp(lengths)/lengths.mean(),np.max(abs(cell.T@cell/np.outer(lengths,lengths)-np.eye(3))))
 ev['primitive_metric_error']=float(metric)
 if metric>.001:ev['reason']='primitive_cell_not_cubic';return None,ev
 H=np.meshgrid(*[np.rint(np.fft.fftfreq(m)*m).astype(int) for m in shape],indexing='ij')
 keep=np.ones(tuple(shape),bool)
 for h,r in zip(H,repeats):keep &= h%r==0
 fractions=[]
 for v in a:
  p=abs(np.fft.fftn(v-v.mean()))**2
  fractions.append(float(p[~keep].sum()/max(p.sum(),1e-30)))
 ev['discarded_power_fraction']=fractions
 if max(fractions)>POLICY['maximum_discarded_power_fraction']:
  ev['reason']='nonperiodic_power_exceeds_limit';return None,ev
 primitive=describe(small,cell);ev['primitive_geometry']=primitive
 if primitive is None or 'SC' not in primitive['geometry_matches'] or primitive['resolved_allowed_shells']['SC']<3:
  ev['reason']='primitive_SC_geometry_or_spectrum_unverified';return None,ev
 # Check unprojected input, including extra domains and unequal domain sizes.
 ci=int(np.argmin([v.mean() for v in a]));v=a[ci]
 template=np.array(list(itertools.product(*[np.arange(r)/r for r in repeats])))
 checks=[]
 for t in POLICY['relative_levels']:
  d=domains(v,t,A,allowed_counts=tuple(range(1,POLICY['maximum_domains']+1)))
  checks.append(d)
  if d['count']!=n or not d.get('isolated') or any(e!=1 for e in d.get('eulers',[])) or d.get('size_cv',1)>.15 or d.get('axis_ratio',1e9)>4:
   ev.update(reason='original_domains_not_equivalent_solid_balls',original_component=ci,original_checks=checks);return None,ev
  c=np.array(d['centers']);best=float('inf')
  for origin in c:
   delta=(c[:,None,:]-origin-template[None,:,:]+.5)%1-.5
   distance=np.linalg.norm(delta@A.T,axis=-1)/lengths.mean()
   i,j=linear_sum_assignment(distance);best=min(best,float(distance[i,j].max()))
  d['SC_repeat_grid_error_primitive_fraction']=best
  if best>POLICY['center_tolerance_primitive_fraction']:
   ev.update(reason='original_center_grid_mismatch',original_component=ci,original_checks=checks);return None,ev
 ev.update(accepted=True,reason='SC_repeat_geometry_and_small_nonperiodic_power_supported',original_component=ci,original_checks=checks,repeats=repeats.tolist(),
  interpretation='approximate SC morphology supported; small nonperiodic component retained as a diagnostic')
 return (small,cell),ev
