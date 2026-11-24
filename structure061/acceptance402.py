"""Conservative checks for a proposed phase, never a phase-name override.

Only retracts an otherwise accepted morphology. No sample IDs, seed labels,
paths, energies or database records enter these checks. Numerical tolerances
are screening policies, not error bounds or thermodynamic stability proofs.
"""
import itertools
import numpy as np
from scipy.optimize import linear_sum_assignment
from particle402 import domains, POLICY as PARTICLE_POLICY
from structure_checks import topology, centers, neighbors

POLICY = dict(relative_levels=[.35, .5, .65], center_tolerance=.02,
              centering_power_tolerance=1e-8, cubic_power_tolerance=1e-3,
              metric_tolerance=1e-3, maximum_size_cv=.15,
              maximum_axis_ratio=4., maximum_domains=64)
FCC = np.array([[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0.]])
PROTOTYPES = {
 'SC': np.array([[0,0,0.]]), 'BCC': np.array([[0,0,0],[.5,.5,.5]]),
 'FCC': FCC, 'DSc': np.vstack([FCC, (FCC+.25)%1]),
 'A15': np.array([[0,0,0],[.5,.5,.5],[0,.5,.25],[0,.5,.75],
                  [.25,0,.5],[.75,0,.5],[.5,.25,0],[.5,.75,0]]),
 'HCP': np.array([[0,0,0],[.5,.5,0],[0,2/3,.5],[.5,1/6,.5]]),
 'iHPa': np.array([[0,0,0],[.5,.5,0.]]),
 'PtS': np.array([[0,0,0],[0,0,.5],[.5,0,.25],[0,.5,.75]])}

def periodic_euler(mask):
    """Euler characteristic of the periodic union of closed voxel cubes.

    Counts vertices/edges/faces/voxels on the torus, not the cut box. Uses
    cubical (26-neighbor foreground) topology; a passing value is necessary
    evidence only, not an identification of a network's skeleton graph.
    """
    d=mask.ndim
    def join(axes):
        v=mask.copy()
        for axis in axes:v |= np.roll(v,1,axis=axis)
        return int(v.sum())
    return sum((-1)**(d-len(axes))*join(axes) for n in range(d+1)
               for axes in itertools.combinations(range(d),n))

def _metric(A):
    length=np.linalg.norm(A,axis=0)
    return length, A.T@A/np.outer(length,length)

def _compatible(A,name):
    length,cos=_metric(A);tol=POLICY['metric_tolerance']
    if not np.allclose(cos,np.eye(3),atol=tol):return False
    if name in ('SC','BCC','FCC','DSc','A15'):
        return np.ptp(length)/length.mean()<tol
    if name in ('HCP','iHPa'):
        return abs(length[1]/length[0]-np.sqrt(3))<tol
    return name=='PtS' and abs(length[0]/length[1]-1)<tol

def prototype_fit(c,A,name):
    """Translation/sign/axis invariant motif fit in a compatible cell metric.

    Axial ratio is free for HCP/iHPa and tetragonal PtS; no ideal close-packing
    c/a is imposed. A bounded conventional-cell model deliberately leaves
    unverified basis/supercell representations for review.
    """
    T=PROTOTYPES[name];c=np.asarray(c)
    if len(c)!=len(T):return None
    best=float('inf')
    for perm in itertools.permutations(range(3)):
        B=A[:,perm]
        if not _compatible(B,name):continue
        cp=c[:,perm]
        for sign in itertools.product((-1,1),repeat=3):
            points=(cp*np.asarray(sign))%1
            for anchor in T:
                delta=(points-points[0])[:,None,:]-(T-anchor)[None,:,:]
                delta-=np.round(delta)
                dist=np.linalg.norm(delta@B.T,axis=-1)/np.mean(np.linalg.norm(B,axis=0))
                i,j=linear_sum_assignment(dist);best=min(best,float(dist[i,j].max()))
    return best if np.isfinite(best) else None

def cubic_power_error(fields):
    """Necessary cubic point-group condition on full signed low-order powers.

    Full intensity is origin independent. This checks rotations/reflections,
    not just equal radial shell positions. Phase information is checked
    separately by the real-space point motif. Nyquist planes are excluded.
    """
    a=np.asarray(fields[0]);k=min(8,min(a.shape)//2-1)
    if a.ndim!=3 or k<2:return None
    h=np.stack(np.meshgrid(*[np.arange(-k,k+1)]*3,indexing='ij'),axis=-1).reshape(-1,3)
    largest=0.
    for v in fields:
        power=abs(np.fft.fftn(v-v.mean()))**2
        p=power[tuple((h%np.array(v.shape)).T)];den=max(np.linalg.norm(p),1e-30)
        for perm in itertools.permutations(range(3)):
            for sign in itertools.product((-1,1),repeat=3):
                g=h[:,perm]*sign
                q=power[tuple((g%np.array(v.shape)).T)]
                largest=max(largest,float(np.linalg.norm(q-p)/den))
    return largest

def _particle_fields(r,fields,A):
    rep=r.get('particle_geometry_representation',{}).get('selection')
    if rep in ('verified_axis_primitive_cell','approximate_axis_primitive_cell'):
        from supercell402 import reduce_supercell
        kw=dict(power_threshold=1e-8,relative_tolerance=1e-3) if rep.startswith('approximate') else {}
        small,B,info=reduce_supercell(fields,A,**kw)
        if info.get('reduced'):return small,B
    return fields,A

def particle_check(r,fields,A,name):
    fields,A=_particle_fields(r,fields,A);A=np.asarray(A)
    if A.shape!=(3,3):return dict(passed=False,reason='particle_cell_not_3d')
    ev=r.get('particle_evidence') or {};checks=[];matches=[]
    # A sphere-forming component need not be the minority component. At large
    # filling fraction the majority component may still be isolated spheres.
    order=sorted(range(len(fields)),key=lambda i:np.mean(fields[i]))
    for ci in order:
        rows=[]
        for level in POLICY['relative_levels']:
            d=domains(np.asarray(fields[ci]),level,A,allowed_counts=tuple(range(1,65)))
            d['prototype_error']=prototype_fit(d.get('centers',[]),A,name)
            solid=(d.get('isolated') and all(x==1 for x in d.get('eulers',[]))
                   and d.get('axis_ratio',1e9)<POLICY['maximum_axis_ratio'])
            equal=(name=='A15' or d.get('size_cv',1)<POLICY['maximum_size_cv'])
            d['passed']=bool(solid and equal and d['prototype_error'] is not None
                             and d['prototype_error']<POLICY['center_tolerance'])
            rows.append(d)
            if not d['passed']:break
        checks.append(dict(component=ci,rows=rows))
        if len(rows)==3 and all(x['passed'] for x in rows):matches.append(ci);break
    evidence=dict(passed=bool(matches),reason='relative_solid_domains_and_point_motif' if matches else 'relative_particle_motif_not_verified',
                  matching_components=matches,components=checks)
    if name in ('BCC','FCC','DSc'):
        # Match the entire centering selection rule, not only the first ten
        # positive-octant shells. Missing allowed motif peaks are not vetoed.
        h=np.meshgrid(*[np.rint(np.fft.fftfreq(n)*n).astype(int) for n in fields[0].shape],indexing='ij')
        allow=(sum(h)%2==0) if name=='BCC' else (h[0]%2==h[1]%2)&(h[1]%2==h[2]%2)
        fractions=[]
        for v in fields:
            p=abs(np.fft.fftn(v-v.mean()))**2
            fractions.append(float(p[~allow].sum()/max(p.sum(),1e-30)))
        evidence['centering_forbidden_fractions']=fractions
        if max(fractions)>POLICY['centering_power_tolerance']:
            evidence.update(passed=False,reason='full_signed_centering_power_conflict')
    if name in ('SC','BCC','FCC','DSc','A15'):
        error=cubic_power_error(fields);evidence['cubic_power_error']=error
        if error is None or error>POLICY['cubic_power_tolerance']:
            evidence.update(passed=False,reason='cubic_directional_power_conflict')
    return evidence

def network_check(fields,A,name):
    # These signatures are scoped to the conventional cells used by the
    # original reference library. They are not universal across unit cells.
    length,cos=_metric(np.asarray(A))
    if len(length)!=3 or np.ptp(length)/length.mean()>1e-3 or not np.allclose(cos,np.eye(3),atol=1e-3):
        return dict(passed=False,reason='network_conventional_cell_not_verified')
    rows=[]
    for ci,v in enumerate(fields):
        checks=[]
        for level in POLICY['relative_levels']:
            m=v>v.min()+level*np.ptp(v);_,_,t=topology(m)
            checks.append(dict(level=level,euler=periodic_euler(m),**t))
        rows.append(dict(component=ci,rows=checks))
    expected_euler={'SG':-4,'DG':-8,'DD':-2}[name]
    counts=[tuple(x['domains'] for x in r['rows']) for r in rows]
    correct_counts=(all(c==(1,1,1) for c in counts) if name!='DG'
                    else sorted(counts)==[(1,1,1),(2,2,2)])
    passed=(len(rows)==2 and correct_counts and all(x['euler']==expected_euler and
            x['spanning_axes']==[0,1,2] for r in rows for x in r['rows']))
    return dict(passed=bool(passed),reason='network_components_and_periodic_euler_supported' if passed else 'network_topology_not_persistent',
                components=rows,expected_periodic_solid_euler=expected_euler,
                scope='necessary reference-cell topology; not a sufficient skeleton identity proof')

def planar_check(fields,A,name):
    from geometry061 import planar
    if len(A)!=2:return dict(passed=False,reason='planar_cell_not_verified')
    if name in ('HEX','C3'):
        e=planar(fields,A,relative=True)
        passed=bool(e and {x['candidate'] for x in e}=={name})
        return dict(passed=passed,reason='relative_planar_domains_supported' if passed else 'relative_planar_geometry_conflict',components=e)
    length,cos=_metric(np.asarray(A));square=(np.ptp(length)/length.mean()<1e-3 and np.allclose(cos,np.eye(2),atol=1e-3))
    checks=[]
    for ci,v in enumerate(fields):
        rows=[]
        for t in POLICY['relative_levels']:
            labs,ids,top=topology(v>v.min()+t*np.ptp(v));c=centers(labs,ids)
            nn=neighbors(c,A) if len(c) else []
            rows.append(dict(level=t,**top,neighbors=nn))
        checks.append(dict(component=ci,rows=rows))
    passed=square and any(all(t['isolated'] and t['domains']==1 and t['neighbors']==[4] for t in r['rows']) for r in checks)
    return dict(passed=bool(passed),reason='primitive_square_domains_supported' if passed else 'square_motif_or_representation_unverified',components=checks)

def apply(result,fields,A):
    r=dict(result);name=r.get('candidate');status=r.get('status','')
    if not name or status.startswith('review_') or name=='DIS':return r
    name={'LAMs':'LAM','PtS1':'PtS','PtS2':'PtS'}.get(name,name)
    if name=='LAM':
        d=r.get('lam_diagnostic')
        evidence=dict(passed=not(d and d.get('family_supported') is False),
                      reason='layer_diagnostic_conflict' if d and d.get('family_supported') is False else 'existing_layer_and_repeat_evidence_retained')
    elif name in PROTOTYPES:evidence=particle_check(r,fields,A,name)
    elif name in ('SG','DG','DD'):evidence=network_check(fields,A,name)
    elif name in ('HEX','C3','SPC'):evidence=planar_check(fields,A,name)
    else:
        evidence=dict(passed=False,reason='dedicated_variant_or_morphology_identity_not_verified')
    r['acceptance_evidence']=dict(policy=POLICY,phase_tested=name,**evidence)
    if not evidence['passed']:
        r.update(candidate=None,status='review_acceptance_'+evidence['reason'],
                 withheld_morphology_candidate=name,withheld_morphology_status=status)
    return r
