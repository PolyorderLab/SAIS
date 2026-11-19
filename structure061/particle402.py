"""Evidence-constrained cubic particle matching.

Allowed reflections may have motif-dependent zeros. Such missing peaks are
waived only with independent centering, isolated genus-zero domains and a
multi-shell spectrum. Forbidden power is NOT waived by this rule.
No input names, sample identifiers, seed labels or database labels are used.
"""
import itertools
import numpy as np
from scipy.ndimage import label
from scipy.optimize import linear_sum_assignment
from skimage.measure import euler_number

TEMPLATES={'SC':np.array([[0,0,0.]]),
 'BCC':np.array([[0,0,0],[.5,.5,.5]]),
 'FCC':np.array([[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0]])}
POLICY=dict(relative_levels=[.35,.5,.65],metric_tolerance=.001,
 center_tolerance=.02,centering_power_tolerance=1e-8,minimum_shells=3,
 maximum_size_cv=.15,maximum_axis_ratio=4.)

def domains(v,level,A,allowed_counts=(1,2,4)):
 mask=v>v.min()+level*np.ptp(v)
 labs,n=label(mask);parent=np.arange(n+1)
 def root(x):
  while parent[x]!=x:
   parent[x]=parent[parent[x]];x=parent[x]
  return x
 for axis in range(3):
  a=np.take(labs,0,axis).ravel();b=np.take(labs,-1,axis).ravel()
  for x,y in set(zip(a[(a>0)&(b>0)],b[(a>0)&(b>0)])):
   parent[root(x)]=root(y)
 labs=np.array([root(i) for i in range(n+1)])[labs]
 ids=np.unique(labs);ids=ids[ids>0]
 result=dict(level=level,count=len(ids),fits={})
 if len(ids) not in allowed_counts:return result
 centers=[];sizes=[];eulers=[];axis_ratios=[];isolated=True
 for ident in ids:
  m=labs==ident;coords=np.argwhere(m)/np.array(v.shape)
  center=(np.angle(np.exp(2j*np.pi*coords).mean(axis=0))/(2*np.pi))%1
  centered=np.roll(m,tuple(np.array(v.shape)//2-np.rint(center*np.array(v.shape)).astype(int)),axis=(0,1,2))
  isolated &= all(not np.take(centered,j,axis).any() for axis in range(3) for j in (0,-1))
  eulers.append(int(euler_number(centered,connectivity=1)))
  xyz=((coords-center+.5)%1-.5)@A.T
  eig=np.linalg.eigvalsh(xyz.T@xyz/max(len(xyz),1))
  axis_ratios.append(float(eig[-1]/max(eig[0],1e-15)))
  centers.append(center);sizes.append(len(coords))
 centers=np.array(centers)
 for name,T in TEMPLATES.items():
  if len(T)!=len(centers):continue
  best=1.
  for origin in centers:
   delta=(centers[:,None,:]-origin-T[None,:,:]+.5)%1-.5
   distance=np.linalg.norm(delta@A.T,axis=-1)/np.mean(np.linalg.norm(A,axis=0))
   i,j=linear_sum_assignment(distance);best=min(best,float(distance[i,j].max()))
  result['fits'][name]=best
 result.update(centers=centers.tolist(),isolated=bool(isolated),eulers=eulers,
  size_cv=float(np.std(sizes)/np.mean(sizes)),axis_ratio=max(axis_ratios))
 return result

def describe(fields,A):
 A=np.asarray(A);a=[np.asarray(x) for x in fields]
 if A.shape!=(3,3) or not a or any(x.ndim!=3 or not np.isfinite(x).all() for x in a):return None
 lengths=np.linalg.norm(A,axis=0)
 metric=max(np.ptp(lengths)/lengths.mean(),np.max(abs(A.T@A/np.outer(lengths,lengths)-np.eye(3))))
 if metric>POLICY['metric_tolerance']:return None
 # Two-component incompressible systems have a unique minority-domain field.
 if len(a)!=2:return None
 v=a[int(np.argmin([x.mean() for x in a]))]
 if np.ptp(v)<1e-5:return None
 checks=[domains(v,t,A) for t in POLICY['relative_levels']]
 matches=[name for name in TEMPLATES if all(
  t['fits'].get(name,1)<POLICY['center_tolerance'] and
  t.get('isolated',False) and all(e==1 for e in t.get('eulers',[])) and
  t.get('size_cv',1)<POLICY['maximum_size_cv'] and
  t.get('axis_ratio',1e9)<POLICY['maximum_axis_ratio'] for t in checks)]
 # Center arrangement is also retained when morphology is distorted.
 arrangements=[name for name in TEMPLATES if all(t['fits'].get(name,1)<POLICY['center_tolerance'] for t in checks)]
 h=np.meshgrid(*[np.rint(np.fft.fftfreq(n)*n).astype(int) for n in v.shape],indexing='ij')
 allow={'BCC':sum(h)%2==0,'FCC':(h[0]%2==h[1]%2)&(h[1]%2==h[2]%2),'SC':np.ones(v.shape,bool)}
 fractions={};shellcounts={}
 q2=sum(x*x for x in h)
 for name,mask in allow.items():
  fractions[name]=max(float((abs(np.fft.fftn(x-x.mean()))**2)[~mask].sum()/max((abs(np.fft.fftn(x-x.mean()))**2).sum(),1e-30)) for x in a)
  p=abs(np.fft.fftn(v-v.mean()))**2
  shells=np.bincount(q2.ravel(),weights=np.where(mask,p,0).ravel());shells[0]=0
  shellcounts[name]=int(np.count_nonzero(shells>max(shells.max(),1e-30)*1e-4))
 return dict(policy=POLICY,metric_error=float(metric),relative_topology=checks,
  arrangements=arrangements,geometry_matches=matches,forbidden_power=fractions,
  resolved_allowed_shells=shellcounts)


def network_evidence(fields):
 from structure_checks import topology
 checks=[]
 for i,v in enumerate(fields):
  rows=[]
  for level in POLICY['relative_levels']:
   _,_,t=topology(v>v.min()+level*np.ptp(v))
   rows.append(dict(level=level,**t))
  checks.append(dict(component=i,checks=rows))
 # Necessary bicontinuity test only; scattering still determines SG symmetry.
 passed=len(checks)>=2 and all(all(t['domains']==1 and t['spanning_axes']==[0,1,2] for t in c['checks']) for c in checks)
 return dict(passed=passed,components=checks,scope='necessary periodic bicontinuity; not sufficient for gyroid identity')

def apply_in_cell(result,fields,A):
 r=dict(result);candidate=r.get('candidate');status=r.get('status','')
 sg=(candidate=='SG' or r.get('scattering_candidate')=='SG')
 network=network_evidence(fields) if sg else None
 if network is not None:r['SG_network_evidence']=network
 ev=describe(fields,A)
 if ev is not None:r['particle_evidence']=ev
 if candidate=='SC' and (ev is None or 'SC' not in ev['geometry_matches']):
  r.update(candidate=None,status='review_SC_center_arrangement_unverified',
   particle_decision='SC_requires_single_equivalent_isolated_solid_ball_per_cubic_cell_at_all_relative_levels')
  return r
 if candidate=='SG' and not network['passed']:
  r.update(candidate=None,status='review_SG_network_topology_conflict',
   particle_decision='SG_requires_periodic_bicontinuous_network_not_disconnected_balls')
 if ev is None:return r
 supported=[n for n in ev['geometry_matches'] if n in ('BCC','FCC') and
  ev['forbidden_power'][n]<=POLICY['centering_power_tolerance'] and
  ev['resolved_allowed_shells'][n]>=POLICY['minimum_shells']]
 if len(supported)!=1:return r
 name=supported[0]
 eligible=(candidate in ('BCC','FCC',None) and
  status in ('review_reflection_or_threshold_conflict','review_scattering_ambiguity',
   'review_particle_candidate_but_nonisolated_domains','review_insufficient_peaks'))
 # A rejected network hypothesis may be resolved by independent sphere evidence.
 network_conflict=sg and not network['passed'] and (candidate=='SG' or status.startswith('review_'))
 if eligible or network_conflict:
  r.update(candidate=name,status='particle_centering_and_topology_supported',
   particle_decision='weak_allowed_reflections_checked_against_centering_and_relative_topology',
   prior_scattering_candidate=candidate,prior_scattering_status=status)
 return r


def apply(result,fields,A):
 # Scattering may keep the conventional cell on a tied score. Geometry must
 # nevertheless be allowed to use an independently verified repeat reduction.
 original=describe(fields,A)
 reduced,cell,validation=fields,A,{}
 if original is None or not original['geometry_matches']:
  from supercell402 import reduce_supercell
  reduced,cell,validation=reduce_supercell(fields,A)
  if validation.get('reduced'):
   primitive=describe(reduced,cell)
   if primitive is not None and len(primitive['geometry_matches'])==1:
    r=apply_in_cell(result,reduced,cell)
    r['particle_geometry_representation']=dict(selection='verified_axis_primitive_cell',
     validation=validation,original_geometry=original,
     scope='geometry only; scattering and physical quality gates retain their input evidence')
    return r
 r=apply_in_cell(result,fields,A)
 if result.get('candidate')=='SC' and not result.get('status','').startswith('review_') and r.get('candidate') is None:
  from approximate_particle402 import evaluate
  proposal,approx=evaluate(fields,A)
  r['approximate_particle_supercell']=approx
  if proposal is not None:
   small,primitive_cell=proposal
   r['candidate']='SC';r['status']='approximate_SC_supercell_geometry_supported'
   r['particle_decision']='SC_supported_by_original_repeat_grid_and_bounded_nonperiodic_power'
   r['particle_evidence']=approx['primitive_geometry']
   r['particle_geometry_representation']=dict(selection='approximate_axis_primitive_cell',validation=approx)
 if validation and not r.get("approximate_particle_supercell",{}).get("accepted"):
  r['particle_geometry_representation']=dict(selection='saved_cell',
   validation=validation,reason='no_unique_verified_primitive_particle_geometry')
 return r
