from run import *
from scipy.signal import find_peaks
def describe(fields,A):
    output=[]
    for ci,v in enumerate(fields):
        if v.ndim!=3:continue
        checks=[]
        for level in (.35,.5,.65):
            labs,ids,t=topology(v>level);row=dict(level=level,**t)
            if len(t['spanning_axes'])==2:
                axis=next(i for i in range(3) if i not in t['spanning_axes']);keep=t['spanning_axes'];profile=v.mean(axis=tuple(keep));n=len(profile)
                peaks,_=find_peaks(np.tile(profile,3),prominence=.05*np.ptp(profile));peaks=peaks[(peaks>=n)&(peaks<2*n)]-n
                B=np.linalg.cholesky(A[:,keep].T@A[:,keep]).T;planes=[]
                for z in peaks:
                    plane=np.take(v,int(z),axis=axis);hl,hi,ht=topology(plane<level)
                    planes.append(dict(index=int(z),holes=len(hi),isolated_holes=ht['isolated'],hole_neighbors=neighbors(centers(hl,hi),B) if ht['isolated'] else []))
                registry=[]
                for z,w in zip(peaks,np.roll(peaks,-1)):
                    a=np.take(v,int(z),axis=axis);b=np.take(v,int(w),axis=axis);a=a-a.mean();b=b-b.mean();den=np.linalg.norm(a)*np.linalg.norm(b)
                    corr=np.fft.ifft2(np.fft.fft2(a)*np.conj(np.fft.fft2(b))).real/max(den,1e-30);ix=np.unravel_index(corr.argmax(),corr.shape)
                    registry.append(dict(from_index=int(z),to_index=int(w),zero_shift=float(corr[0,0]),best=float(corr[ix]),shift_grid=list(map(int,ix)),shape=list(a.shape)))
                row.update(normal_axis=axis,planes=planes,registry=registry)
            checks.append(row)
        output.append(dict(component='phi'+chr(65+ci),checks=checks))
    return output
if __name__=='__main__':
    out=base/'details';out.mkdir(exist_ok=True)
    cases=json.loads((root/'phase_identifier_structure042/all_cases.json').read_text())['cases']
    for key in ['ref_PLAAA','ref_PLABA','ref_Torus','ref_AB6']:
        f,A,h=load(Path(next(c['input'] for c in cases if c['id']==key)));r=describe(f,A);(out/(key+'.json')).write_text(json.dumps(r,indent=2));print(key,r[0],flush=True)
