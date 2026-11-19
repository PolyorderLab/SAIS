"""Read-only peak/spacing and sampling diagnostics; no classification override."""
import numpy as np
from scipy.signal import find_peaks,peak_widths

def _cv(values):
    a=np.asarray(values,dtype=float)
    return None if len(a)<2 else float(np.std(a)/max(abs(np.mean(a)),1e-15))
def inspect_peaks(fields,A):
    component=int(np.argmin([np.mean(x) for x in fields]))
    v=np.asarray(fields[component],dtype=float)
    if v.ndim!=1 or len(v)<12 or not np.isfinite(v).all() or np.ptp(v)<.01:return None
    n=len(v);L=float(np.linalg.norm(np.asarray(A)[:,0]));dx=L/n;tile=np.tile(v,3)
    indices,props=find_peaks(tile,prominence=.05*np.ptp(v))
    take=(indices>=n)&(indices<2*n);indices=indices[take];prom=props['prominences'][take]
    if not len(indices):return dict(component=component,peak_count=0,status='no_resolved_peaks')
    half,_,left,right=peak_widths(tile,indices,rel_height=.5)
    _,_,l10,r10=peak_widths(tile,indices,rel_height=.1)
    _,_,l90,r90=peak_widths(tile,indices,rel_height=.9)
    spacings=(np.roll(indices,-1)-indices)%n
    if len(indices)==1:spacings[:]=n
    maxprom=float(max(prom));peaks=[]
    for j,i in enumerate(indices):
        peaks.append(dict(position=float((i-n)*dx),height=float(tile[i]),prominence=float(prom[j]),
          relative_prominence=float(prom[j]/maxprom),weak=bool(prom[j]<.5*maxprom),
          half_height_width=float(half[j]*dx),half_height_grid_points=float(half[j]),
          left_interface_grid_points=float(l10[j]-l90[j]),right_interface_grid_points=float(r90[j]-r10[j]),
          next_peak_distance=float(spacings[j]*dx)))
    minimum_width=float(min(half));min_interface=float(min(np.min(l10-l90),np.min(r90-r10)))
    sampling='coarse' if minimum_width<4 or min_interface<2 else 'borderline' if minimum_width<6 or min_interface<3 else 'resolved_by_screen'
    pcv=_cv(prom);scv=_cv(spacings);wcv=_cv(half)
    differences=[]
    if any(p['weak'] for p in peaks):differences.append('strong_weak_peaks')
    if scv is not None and scv>.05:differences.append('unequal_spacing')
    if pcv is not None and pcv>.05:differences.append('unequal_prominence')
    if wcv is not None and wcv>.05:differences.append('unequal_width')
    return dict(component=component,component_name='phiA' if component==0 else 'phiB',grid_points=n,grid_spacing=dx,
      peak_count=len(peaks),weak_peak_count=sum(p['weak'] for p in peaks),peaks=peaks,
      prominence_cv=pcv,spacing_cv=scv,width_cv=wcv,differences=differences,
      sampling=sampling,min_half_width_grid_points=minimum_width,min_interface_grid_points=min_interface,
      rules=dict(min_prominence_fraction=.05,weak_relative_prominence=.5,variation_notice_cv=.05,
                 coarse_half_width_points=4,coarse_interface_points=2),
      interpretation='Screening diagnostics only; interpolation does not establish grid convergence.')
