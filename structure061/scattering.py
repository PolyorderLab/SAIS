"""Experimental Fourier-only reference matching; no sample labels or geometry models."""
import itertools
import numpy as np
import h5py
def align_peaks(q,p,expected,tolerance=.01):
    """Monotone alignment: weak insertions have small cost, never shift all peaks.

    Missing expected reflections retain a fixed cost. This is an engineering
    extension, not Eq6. Intensities are normalized within a common q window.
    """
    q=np.asarray(q);p=np.asarray(p);expected=np.asarray(expected)
    m,n=len(q),len(expected)
    weights=p/(p.sum() or 1)
    dp=np.full((m+1,n+1),np.inf);dp[0,0]=0
    for i in range(m+1):
        for j in range(n+1):
            if i<m:dp[i+1,j]=min(dp[i+1,j],dp[i,j]+weights[i])
            if j<n:dp[i,j+1]=min(dp[i,j+1],dp[i,j]+1/max(n,1))
            if i<m and j<n:
                delta=abs(q[i]/expected[j]-1)
                # Existing shell resolution is 0.001: do not rank rounded
                # reference matrices by meaningless sub-resolution differences.
                if delta<=tolerance:dp[i+1,j+1]=min(dp[i+1,j+1],dp[i,j]+weights[i]*(max(0,delta-.001)/tolerance)**2)
    return float(dp[m,n])
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

def references(path):
    def read(f,r):
        a=f[r][...]
        return [read(f,v) for v in a.flat] if a.dtype.kind=='O' else a.T.squeeze()
    out={}
    with h5py.File(path) as f:
        for name in f['dict']:
            pair=f['dict'][name][...].flat; p,m=list(pair)
            rows=[]
            for flag,hs in f[p][...].T:
                groups=read(f,hs)
                if isinstance(groups,np.ndarray):groups=[groups]
                rows.append((float(np.asarray(read(f,flag))),[np.asarray(h,dtype=int).reshape(-1) for h in groups]))
            out[name]=(rows,f[m][...].T)
    return out

def system(A):
    g=A.T@A;e=np.sqrt(np.diag(g));c=g/np.outer(e,e)
    if np.allclose(c,np.eye(len(e)),atol=1e-3):
        return 'cubic' if len(e)==3 and np.ptp(e)/e.mean()<1e-3 else 'orthogonal'
    if len(e)>=2 and abs(e[0]/e[1]-1)<1e-3 and abs(abs(c[0,1])-.5)<1e-3:
        if len(e)==2 or max(abs(c[0,2]),abs(c[1,2]))<1e-3:return 'hexagonal'
    return 'other'

def identify(field,A,refs):
    power=abs(np.fft.fftn(field-field.mean())/field.size)**2
    peak=float(power.max());dim=field.ndim
    if peak<1e-16:return dict(status='no_resolved_scattering',candidate=None)
    grids=np.meshgrid(*[np.rint(np.fft.fftfreq(n)*n).astype(int) for n in field.shape],indexing='ij')
    indices=np.stack(grids,axis=-1)
    transforms=[np.eye(dim,dtype=int)];primitive=[]
    for threshold in (1e-5,1e-4,1e-3):
        h=indices[power>threshold*peak]
        H=np.array(hermite_normal_form(Matrix(h.T)),dtype=int)
        if H.shape==(dim,dim):primitive.append(H)
    if len(primitive)==3 and all(np.array_equal(primitive[0],h) for h in primitive[1:]):
        H=primitive[0]
        if abs(np.linalg.det(H))>1 and abs(np.linalg.det(H))<=64:transforms.append(H)
    trials=[]
    for H in transforms:
        cell=A@np.linalg.inv(H).T;kind=system(cell)
        for name,(pattern,M) in refs.items():
            if M.shape!=(dim,dim) or system(np.linalg.inv(M).T)!=kind:continue
            if kind=='other':continue
            h=np.array(list(itertools.product(range(12),repeat=dim)))
            h=h[np.any(h,axis=1)];mapped=h@H.T
            valid=np.all(abs(mapped)<=np.array(field.shape)//2,axis=1)
            h=h[valid];mapped=mapped[valid]
            q=np.linalg.norm(h@M.T,axis=1);p=power[tuple((mapped%field.shape).T)]
            order=np.argsort(q);shells=[]
            for idx in order:
                if shells and abs(q[idx]-shells[-1][0])<=.001*shells[-1][0]:shells[-1][1]+=float(p[idx])
                else:shells.append([float(q[idx]),float(p[idx])])
            rq=np.array([np.linalg.norm(M@hs[0][:dim]) for flag,hs in pattern]);flags=np.array([flag>0 for flag,hs in pattern])
            limit=rq[min(9,len(rq)-1)]
            for threshold in (1e-6,1e-5,1e-4):
                observed=np.array([s[0] for s in shells if s[1]>threshold*peak])
                usable=[s for s in shells if s[0]<=limit*(1+1e-6) and s[1]>threshold*peak]
                aligned=align_peaks([s[0] for s in usable],[s[1] for s in usable],rq[flags & (rq<=limit*(1+1e-6))])
                expected=rq[flags];n=min(10,len(observed),len(expected))
                eps=float(np.linalg.norm(1/observed[:n]-1/expected[:n])/n) if n else None
                checks=[]
                for qref,flag in zip(rq[:10],flags[:10]):
                    matches=[s[1] for s in shells if abs(s[0]-qref)<=.001*qref]
                    value=sum(matches);checks.append(dict(q=float(qref),expected=bool(flag),relative_power=value/peak,observed=bool(value>threshold*peak)))
                trials.append(dict(phase=name,transform=H.tolist(),cell_system=kind,threshold=threshold,
                    alignment_score=aligned,epsilon=eps,Nq=n,observed_q=observed[:10].tolist(),expected_q=expected[:10].tolist(),mismatches=sum(x['expected']!=x['observed'] for x in checks),checks=checks))
    # Compare the same number of peaks for every candidate/threshold. A short
    # spectrum is not an ambiguity; it is insufficient evidence below four peaks.
    common_n=min((x['Nq'] for x in trials),default=0)
    for row in trials:
        n=common_n
        row['comparison_Nq']=n
        row['comparison_epsilon']=float(np.linalg.norm(1/np.array(row['observed_q'][:n])-1/np.array(row['expected_q'][:n]))/n) if n else None
    ranked=[]
    for name in refs:
        for H in transforms:
            rows=[x for x in trials if x['phase']==name and x['transform']==H.tolist() and common_n>=4]
            if len(rows)!=3:continue
            checks=rows[0]['checks']
            forbidden=sum(c['relative_power'] for c in checks if not c['expected'])
            total=sum(c['relative_power'] for c in checks)
            ranked.append(dict(phase=name,transform=H.tolist(),comparison_Nq=common_n,
                forbidden_shell_power_fraction=forbidden/total if total else None,
                strongest_forbidden_relative_power=max((c['relative_power'] for c in checks if not c['expected']),default=0),
                median_epsilon=float(np.median([x['comparison_epsilon'] for x in rows])),
                alignment_scores=[x['alignment_score'] for x in rows],median_alignment=float(np.median([x['alignment_score'] for x in rows])),
                mismatches=[x['mismatches'] for x in rows],epsilons=[x['epsilon'] for x in rows]))
    ranked.sort(key=lambda x:(x['median_alignment'],x['median_epsilon'],x['phase']))
    tied=[x for x in ranked if abs(x['median_alignment']-ranked[0]['median_alignment'])<1e-12 and abs(x['median_epsilon']-ranked[0]['median_epsilon'])<1e-12] if ranked else []
    candidate=ranked[0]['phase'] if len({x['phase'] for x in tied})==1 else None
    status=('review_insufficient_peaks' if common_n<4 else 'review_scattering_ambiguity') if candidate is None else 'scattering_supported' if max(ranked[0]['mismatches'])==0 else 'review_reflection_or_threshold_conflict'
    return dict(status=status,candidate=candidate,selection_rule='intensity-weighted monotone peak alignment; Eq6 secondary diagnostic; experimental extension',
        ranking=ranked,primitive_transforms=[h.tolist() for h in transforms],trials=trials)
