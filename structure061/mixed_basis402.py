"""Integer mixed-basis search with all-component reciprocal-lattice validation.
Candidate cells are selected by metric, never by sample name or phase label.
"""
import itertools
import numpy as np
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form

def primitive_basis(fields,A):
    shape=fields[0].shape;d=len(shape)
    hs=np.stack(np.meshgrid(*[np.rint(np.fft.fftfreq(n)*n).astype(int) for n in shape],indexing='ij'),axis=-1)
    active=np.zeros(shape,dtype=bool)
    for x in fields:
        p=abs(np.fft.fftn(x-x.mean())/x.size)**2
        if p.max()>1e-24:active |= p>p.max()*1e-12
    h=hs[active]
    if not len(h):return None,None
    H=np.array(hermite_normal_form(Matrix(h.T)),dtype=int)
    if H.shape!=(d,d):return None,None
    T=np.linalg.inv(H).T
    B=np.asarray(A)@T
    # Integer size reduction: shorten the primitive basis before enumerating.
    for _ in range(60):
        order=np.argsort(np.sum(B*B,axis=0));B=B[:,order];T=T[:,order]
        changed=False
        for j in range(1,d):
            for i in range(j-1,-1,-1):
                q=int(np.rint(B[:,i]@B[:,j]/(B[:,i]@B[:,i])))
                if q:B[:,j]-=q*B[:,i];T[:,j]-=q*T[:,i];changed=True
        if not changed:break
    return B,T

def candidate_cells(fields,A,max_candidates=18):
    from scattering import system
    B,T=primitive_basis(fields,A)
    if B is None:return []
    d=len(A)
    coeff=np.array([v for v in itertools.product(range(-2,3),repeat=d)
                    if any(v) and next(z for z in v if z)!=0 and next(z for z in v if z)>0],dtype=int)
    vectors=coeff@B.T;lens=np.linalg.norm(vectors,axis=1)
    order=np.argsort(lens);coeff=coeff[order];vectors=vectors[order];lens=lens[order]
    keep=lens<=lens[0]*3.01;coeff=coeff[keep];vectors=vectors[keep];lens=lens[keep]
    products=vectors@vectors.T/np.outer(lens,lens);found={}
    for ids in itertools.combinations(range(len(coeff)),d):
        vs=vectors[list(ids)].T;K=coeff[list(ids)].T
        det=abs(round(np.linalg.det(K)))
        if not det or det>16:continue
        # Only primitive/conventional cells with supported metric families.
        if d==2:
            cos=products[ids[0],ids[1]]
            if abs(cos)<1e-3:pass
            elif abs(abs(cos)-.5)<1e-3 and abs(lens[ids[0]]/lens[ids[1]]-1)<1e-3:
                if cos>0:K[:,1]*=-1;vs[:,1]*=-1
            else:continue
        else:
            pairs=[(0,1,2),(0,2,1),(1,2,0)]
            diagonal=all(abs(products[ids[i],ids[j]])<1e-3 for i,j in [(0,1),(0,2),(1,2)])
            if not diagonal:
                match=None
                for i,j,k in pairs:
                    if abs(lens[ids[i]]/lens[ids[j]]-1)<1e-3 and abs(abs(products[ids[i],ids[j]])-.5)<1e-3 and max(abs(products[ids[k],ids[i]]),abs(products[ids[k],ids[j]]))<1e-3:
                        match=(i,j,k);break
                if match is None:continue
                K=K[:,match];vs=vs[:,match]
                if vs[:,0]@vs[:,1]>0:K[:,1]*=-1;vs[:,1]*=-1
        kind=system(vs)
        if kind=='other':continue
        # Rotation/reflection-equivalent metrics share the same trial.
        gram=vs.T@vs
        key=tuple(np.round(gram.ravel(),8))
        if key not in found:found[key]=(det,kind,vs,T@K)
    ordered=sorted(found.values(),key=lambda z:(z[0],0 if z[1] in ('cubic','hexagonal') else 1,tuple(np.diag(z[2].T@z[2]))))
    return ordered[:max_candidates]

def transform_fields(fields,T,relative_tolerance=1e-6,max_voxels=8000000):
    shape=fields[0].shape;d=len(shape)
    h=np.stack(np.meshgrid(*[np.rint(np.fft.fftfreq(n)*n).astype(int) for n in shape],indexing='ij'),axis=-1)
    mapped=h@T;integer=np.rint(mapped).astype(int)
    keep=np.max(abs(mapped-integer),axis=-1)<1e-7
    Fs=[np.fft.fftn(x)/x.size for x in fields]
    errors=[]
    for x,F in zip(fields,Fs):
        reconstructed=np.fft.ifftn(np.where(keep,F,0)).real*x.size
        errors.append(float(np.max(abs(x-reconstructed)))/max(float(np.ptp(x)),1e-12))
    if max(errors)>relative_tolerance:return None,dict(accepted=False,reason='nonperiodic_components',errors=errors)
    # No nonintegral reflection is folded into another index.
    mapped=integer[keep]
    n=2*np.max(abs(mapped),axis=0)+1
    if np.prod(n)>max_voxels:return None,dict(accepted=False,reason='grid_budget',shape=n.tolist())
    targets=tuple((mapped%n).T);out=[]
    for F in Fs:
        reduced=np.zeros(n,dtype=complex);np.add.at(reduced,targets,F[keep])
        wave=np.fft.ifftn(reduced)*np.prod(n)
        imag=float(np.max(abs(wave.imag)))/max(float(np.ptp(wave.real)),1e-12)
        if imag>relative_tolerance:return None,dict(accepted=False,reason='nonreal_spectral_mapping',relative_imaginary=imag)
        out.append(wave.real)
    return out,dict(accepted=True,errors=errors,shape=n.tolist(),basis_transform=T.tolist())

def metric_compatible_references(B,refs,tolerance=1e-3):
    """Require conventional-cell shape agreement, allowing overall scale.
    Restricted to the mixed-basis fallback; existing distortion policies remain.
    Axes are compared in their actual indexing order, not silently permuted.
    """
    G=B.T@B;G=G/np.trace(G);out={}
    for name,(pattern,M) in refs.items():
        C=np.linalg.inv(M).T
        if C.shape!=B.shape:continue
        Q=C.T@C;Q=Q/np.trace(Q)
        if np.max(abs(abs(G)-abs(Q)))<=tolerance:out[name]=(pattern,M)
    return out
