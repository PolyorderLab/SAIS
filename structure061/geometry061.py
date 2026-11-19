"""Phase-independent continuous translation and relative-level geometry tests."""
import numpy as np
from scipy.ndimage import binary_fill_holes,label
from structure_checks import topology,centers,neighbors

def lamellar(fields):
    a=fields[0]
    if len(fields)!=2 or a.ndim not in (1,2) or min(np.ptp(x) for x in fields)<1e-4:return None
    F=np.fft.fftn(a-a.mean());power=abs(F)**2
    axes=[np.rint(np.fft.fftfreq(n)*n).astype(int) for n in a.shape]
    H=np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1)
    h=H[np.unravel_index(np.argmax(power),a.shape)]
    g=np.gcd.reduce(abs(h));direction=h//max(g,1)
    # All significant Fourier vectors must lie on one reciprocal-space line.
    on=np.ones(a.shape,dtype=bool) if a.ndim==1 else H[...,0]*direction[1]==H[...,1]*direction[0]
    leak=max(float((abs(np.fft.fftn(x-x.mean()))**2)[~on].sum()/max((abs(np.fft.fftn(x-x.mean()))**2).sum(),1e-30)) for x in fields)
    if leak>1e-8:return None
    phase=np.exp(2j*np.pi*np.sum(H*h,axis=-1)/float(h@h))
    error=max(float(np.max(abs(x-np.fft.ifftn(np.fft.fftn(x)*phase).real)))/np.ptp(x) for x in fields)
    if error>1e-3:return None
    return dict(model='LAM',reciprocal_mode=h.tolist(),off_axis_power_fraction=leak,continuous_translation_error=error)

def planar(fields,A,relative=False):
    if len(A)!=2:return []
    evidence=[]
    for ci,v in enumerate(fields):
        if np.ptp(v)<1e-4:continue
        checks=[]
        for t in (.35,.5,.65):
            level=float(v.min()+t*np.ptp(v)) if relative else t
            labs,ids,top=topology(v>level)
            if not top['isolated']:break
            c=centers(labs,ids);nn=neighbors(c,A);roundness=[];valid=True
            for ident,center in zip(ids,c):
                # Unwrap on a tiled grid: a valid tilted-cell domain can touch
                # the saved array edge even after circular-centroid recentering.
                original=labs==ident
                tiled=np.tile(original,(3,3));labels,_=label(tiled)
                seed=np.argwhere(original)[0]+np.array(v.shape)
                mask=labels==labels[tuple(seed)]
                if any(np.take(mask,k,axis=ax).any() for ax in (0,1) for k in (0,-1)) or np.any(binary_fill_holes(mask)&~mask):valid=False;break
                pos=np.argwhere(mask)/np.array(v.shape)@A.T
                if len(pos)<6:valid=False;break
                eig=np.linalg.eigvalsh(np.cov(pos.T));ratio=float(eig[-1]/max(eig[0],1e-20));roundness.append(ratio)
                if relative and ratio>1.3:valid=False;break
            if not valid:break
            name='C3' if all(n==3 for n in nn) else 'HEX' if all(n==6 for n in nn) else None
            checks.append(dict(level=level,domains=len(ids),neighbors=nn,candidate=name,centers=c.tolist(),covariance_ratios=roundness,hole_free=True))
        if len(checks)==3 and checks[0]['candidate'] and len({x['candidate'] for x in checks})==1:
            evidence.append(dict(component='phi'+chr(65+ci),candidate=checks[0]['candidate'],checks=checks,relative_levels=relative))
    return evidence
