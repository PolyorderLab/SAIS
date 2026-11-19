"""LAM family diagnostics for resolved binary 1D fields.
Morphology, repeat equivalence and thermodynamic stability are separate facts.
"""
import math
import numpy as np

def _segments(v,level,L):
    nxt=np.roll(v,-1);up=[];down=[]
    for i,(a,b) in enumerate(zip(v,nxt)):
        if a<=level<b:up.append((i+(level-a)/(b-a))/len(v))
        elif a>level>=b:down.append((i+(level-a)/(b-a))/len(v))
    if not up or len(up)!=len(down):return None
    up=np.array(up);down=np.array(down)
    aw=np.array([np.min((down-u)%1) for u in up])
    starts=np.sort(up);spacings=(np.roll(starts,-1)-starts)%1
    if len(starts)==1:spacings[:]=1
    bw=spacings-aw
    def cv(a):return None if len(a)<2 else float(np.std(a)/max(np.mean(a),1e-15))
    return dict(layer_count=len(up),A_widths=(aw*L).tolist(),B_widths=(bw*L).tolist(),
                repeat_spacings=(spacings*L).tolist(),A_width_cv=cv(aw),B_width_cv=cv(bw),
                spacing_cv=cv(spacings))
def diagnose(fields,A):
    if len(fields)!=2 or fields[0].ndim!=1:return None
    f=[np.asarray(x,dtype=float) for x in fields];a=f[0]
    if len(a)<12 or any(x.shape!=a.shape or not np.isfinite(x).all() for x in f):return None
    spans=[float(np.ptp(x)) for x in f]
    if min(spans)<.01 or max(float(np.std(x)) for x in f)<.003:return None
    if max(abs(f[0]+f[1]-1))>1e-3 or any(x.min() < -1e-3 or x.max()>1.001 for x in f):return None
    L=float(np.linalg.norm(np.asarray(A)[:,0]))
    masks=[_segments(a,float(a.min()+t*np.ptp(a)),L) for t in [.35,.5,.65]]
    if any(m is None for m in masks):return None
    layers=[m['layer_count'] for m in masks]
    if len(set(layers))!=1 or len(a)/layers[0]<6:return None
    F=[np.fft.fft(x) for x in f];h=np.rint(np.fft.fftfreq(len(a))*len(a)).astype(int)
    power=abs(F[0])**2;power[0]=0
    high=float(power[abs(h)>.4*len(a)].sum()/max(power.sum(),1e-30))
    if high>.01:return None
    dominant=int(abs(h[np.argmax(power)]))
    positive=np.fft.rfft(a-a.mean());p=abs(positive)**2
    prominent=[int(i) for i in np.argsort(p[1:])[::-1][:8]+1 if p[i]>p.max()*1e-6]
    autocorr=np.fft.ifft(abs(F[0]-np.where(h==0,F[0],0))**2).real
    autocorr/=max(float(autocorr[0]),1e-30)
    # All resolvable integer repeat counts are tested. r=1 is deliberately
    # excluded: periodic boundaries alone are not evidence for a short period.
    trials=[]
    for repeats in range(2,min(64,len(a)//6)+1):
        phase=np.exp(2j*np.pi*h/repeats)
        errors=[float(np.max(abs(x-np.fft.ifft(z*phase).real)))/s for x,z,s in zip(f,F,spans)]
        trials.append(dict(repeats=repeats,period=L/repeats,grid_points=len(a)/repeats,
                           normalized_error=max(errors),passed=max(errors)<=1e-3))
    passed=[t for t in trials if t['passed']]
    chosen=max(passed,key=lambda t:t['repeats']) if passed else None
    return dict(morphology_family='LAM',evidence_scope='binary one-dimensional layered density; simulation may constrain dimension',
                dominant_mode=dominant,significant_modes=prominent,cell_length=L,
                layer_checks=masks,high_frequency_power_fraction=high,
                repeat_trials=trials,verified_repeats=chosen['repeats'] if chosen else None,
                verified_period=chosen['period'] if chosen else None,
                period_status='nontrivial_repeat_verified' if chosen else 'short_repeat_unverified',
                whole_cell_shift_not_used_as_validation=True,thermodynamic_stability='not_assessed')

def describe(fields,A):
    result=diagnose(fields,A)
    if result is not None:
        from lam_peak402 import inspect_peaks
        result['peak_diagnostic']=inspect_peaks(fields,A)
        result.update(family_supported=True,repeat_error_tolerance=1e-3,
                      relative_interface_levels=[.35,.5,.65])
        return result
    if len(fields)!=2 or fields[0].ndim!=1:return None
    a=np.asarray(fields[0]);reason='family_quality_checks_not_passed';counts=[]
    if np.isfinite(a).all() and len(a)>=12 and np.ptp(a)>=.01:
        masks=[_segments(a,float(a.min()+t*np.ptp(a)),float(np.linalg.norm(np.asarray(A)[:,0]))) for t in [.35,.5,.65]]
        counts=[m['layer_count'] if m else None for m in masks]
        if None in counts:reason='missing_consistent_interfaces'
        elif len(set(counts))!=1:reason='layer_count_changes_with_threshold'
        elif len(a)/counts[0]<6:reason='insufficient_grid_points_per_layer'
        else:
            F=abs(np.fft.fft(a-a.mean()))**2
            h=np.fft.fftfreq(len(a))*len(a)
            if F[abs(h)>.4*len(a)].sum()/max(F.sum(),1e-30)>.01:reason='high_frequency_content'
    return dict(family_supported=False,morphology_family=None,reason=reason,
                layer_counts=counts,thermodynamic_stability='not_assessed')
