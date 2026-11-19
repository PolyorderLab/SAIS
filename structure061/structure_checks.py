"""Phase-name-independent descriptors and explicit geometric model tests."""
import itertools
import numpy as np
from scipy.ndimage import label,generate_binary_structure
from scipy.signal import find_peaks
LEVELS=(.35,.50,.65)
def periodic_labels(mask):
    labs,n=label(mask,generate_binary_structure(mask.ndim,1));parent=np.arange(n+1)
    def find(a):
        while parent[a]!=a:parent[a]=parent[parent[a]];a=parent[a]
        return a
    for axis in range(mask.ndim):
        a=np.take(labs,0,axis=axis).ravel();b=np.take(labs,-1,axis=axis).ravel();v=(a>0)&(b>0)
        for x,y in np.unique(np.stack([a[v],b[v]],axis=1),axis=0):
            i,j=find(x),find(y)
            if i!=j:parent[j]=i
    res=np.array([find(i) for i in range(n+1)])[labs]
    return res,np.unique(res[res>0])
def topology(mask):
    labs,ids=periodic_labels(mask);n=len(ids)
    doubled=[len(periodic_labels(np.concatenate([mask,mask],axis=i))[1]) for i in range(mask.ndim)]
    return labs,ids,dict(domains=n,doubled_counts=doubled,
        isolated=bool(n and all(x==2*n for x in doubled)),
        spanning_axes=[i for i,x in enumerate(doubled) if n and x==n])
def centers(labs,ids):
    return np.array([np.angle(np.mean(np.exp(2j*np.pi*np.argwhere(labs==i)/np.array(labs.shape)),axis=0))/(2*np.pi)%1 for i in ids])
def neighbors(c,A):
    tiled=(c[None,:,:]+np.array(list(itertools.product([-1,0,1],repeat=A.shape[0])))[:,None,:]).reshape(-1,A.shape[0]);out=[]
    for x in c:
        d=np.linalg.norm((tiled-x)@A.T,axis=1);d=np.sort(d[d>1e-6]);out.append(int(np.sum(d<=d[0]*1.05)))
    return out
def fcc_model(labs,ids,A,field):
    if len(ids)!=4:return None
    G=A.T@A;edges=np.sqrt(np.diag(G));cos=G/np.outer(edges,edges)
    if not np.allclose(cos,np.eye(3),atol=1e-4) or np.ptp(edges)/edges.mean()>1e-3:return None
    c=centers(labs,ids);ref=np.array([[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0]])
    # Translation invariant set comparison in the saved conventional cubic cell.
    errors=[]
    for origin in c:
        shifted=(c-origin)%1;delta=shifted[:,None,:]-ref[None,:,:];delta-=np.round(delta)
        dist=np.linalg.norm(delta@A.T,axis=2)
        errors.append(max(dist.min(axis=0).max(),dist.min(axis=1).max())/edges.mean())
    error=float(min(errors));nn=neighbors(c,A)
    if error>.03 or nn!=[12]*4:return None
    volumes=np.array([np.sum(labs==i) for i in ids],dtype=float);volume_cv=float(volumes.std()/volumes.mean())
    if volume_cv>.05:return None
    power=np.abs(np.fft.fftn(field-field.mean())/field.size)**2
    signed=[np.rint(np.fft.fftfreq(n)*n).astype(int) for n in field.shape];H=np.meshgrid(*signed,indexing='ij')
    parity=[h%2 for h in H];allowed=(parity[0]==parity[1])&(parity[1]==parity[2])
    low=np.logical_and.reduce([abs(h)<=4 for h in H]);total=float(power[low].sum())
    leakage=float(power[low&~allowed].sum()/total) if total>0 else 1.
    if leakage>1e-4:return None
    from skimage.measure import euler_number
    balls=[]
    for i,center in zip(ids,c):
        mask=labs==i;shift=np.array(mask.shape)//2-np.rint(center*np.array(mask.shape)).astype(int)
        mask=np.roll(mask,tuple(shift),axis=(0,1,2));pos=np.argwhere(mask)
        lo=pos.min(axis=0);hi=pos.max(axis=0)+1;crop=mask[tuple(slice(a,b) for a,b in zip(lo,hi))]
        euler=int(euler_number(crop,connectivity=1));balls.append(euler==1)
    if not all(balls):return None
    return dict(model='FCC',centers_fractional=c.tolist(),position_error_cell_fraction=error,neighbors=nn,solid_ball_euler=balls,
                domain_volume_cv=volume_cv,normalized_low_order_forbidden_energy_fraction=leakage,forbidden_energy_limit=1e-4)
def hpl_model(field,labs,ids,top,A,level):
    axes=top['spanning_axes']
    if len(axes)!=2 or not ids.size:return None
    normal=next(i for i in range(3) if i not in axes)
    if top['doubled_counts'][normal]!=2*len(ids):return None
    G=A.T@A;cos=G/np.sqrt(np.outer(np.diag(G),np.diag(G)))
    if max(abs(cos[normal,j]) for j in axes)>1e-4:return None
    plane_basis=A[:,axes];A2=np.linalg.cholesky(plane_basis.T@plane_basis).T;holes=[]
    for i in ids:
        occupancy=np.sum(labs==i,axis=tuple(axes));index=int(np.argmax(occupancy))
        mask=np.take(field,index,axis=normal)<level;hl,hi,ht=topology(mask)
        if not ht['isolated'] or not len(hi):return None
        c=centers(hl,hi);nn=neighbors(c,A2)
        if any(n!=6 for n in nn):return None
        holes.append(dict(slice_index=index,hole_count=len(hi),neighbors=nn))
    return dict(model='HPL',layer_normal_axis=normal,layers=len(ids),hole_planes=holes)
def diamond_sphere_model(labs,ids,A):
    if len(ids)!=8:return None
    G=A.T@A;edges=np.sqrt(np.diag(G))
    if not np.allclose(G/np.outer(edges,edges),np.eye(3),atol=1e-4) or np.ptp(edges)/edges.mean()>1e-3:return None
    c=centers(labs,ids)
    if neighbors(c,A)!=[4]*8:return None
    fcc=np.array([[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0]])
    errors=[]
    from scipy.optimize import linear_sum_assignment
    for sign in (-1,1):
        ref=np.vstack([fcc,(fcc+sign*.25)%1])
        for origin in c:
            delta=(c-origin)[:,None,:]-ref[None,:,:];delta-=np.round(delta)
            dist=np.linalg.norm(delta@A.T,axis=2)/edges.mean()
            i,j=linear_sum_assignment(dist);errors.append(float(dist[i,j].max()))
    if min(errors)>.03:return None
    from skimage.measure import euler_number
    for i,center in zip(ids,c):
        mask=np.roll(labs==i,tuple(np.array(labs.shape)//2-np.rint(center*np.array(labs.shape)).astype(int)),axis=(0,1,2))
        if euler_number(mask,connectivity=1)!=1:return None
    volumes=np.array([(labs==i).sum() for i in ids]);cv=float(volumes.std()/volumes.mean())
    if cv>.05:return None
    return dict(model='DSc',centers_fractional=c.tolist(),neighbors=[4]*8,position_error_cell_fraction=min(errors),domain_volume_cv=cv)
def lam_1d(fields):
    if len(fields)!=2 or fields[0].ndim!=1:return None
    a=fields[0];span=float(np.ptp(a));n=len(a)
    if n<8 or span<1e-4:return None
    for x in fields:
        jumps=np.abs(np.roll(x,-1)-x);typical=float(np.quantile(jumps,.95))
        if jumps.max()>4*max(typical,1e-12):return None
    # Count prominent peaks on three repeated periods to avoid edge-origin bias.
    peaks,_=find_peaks(np.tile(a,3),prominence=.1*span,distance=2);peaks=peaks[(peaks>=n)&(peaks<2*n)]-n
    if not len(peaks):return None
    spacing=np.diff(np.r_[peaks,peaks[0]+n]);period=n/len(peaks)
    if np.max(abs(spacing-period))>1 or abs(period-round(period))>1e-9:return None
    period=int(round(period));error=max(float(np.max(abs(x-np.roll(x,period))))/max(float(np.ptp(x)),1e-12) for x in fields)
    if error>1e-3:return None
    return dict(model='LAM',period_grid_points=period,periods_in_saved_cell=len(peaks),normalized_translation_error=error)
def analyze(fields,A):
    lam=lam_1d(fields)
    if lam:return dict(candidate='LAM',evidence=lam,reason='periodic_nonuniform_one_dimensional_binary_density')
    if fields[0].ndim!=3:return dict(candidate=None,reason='no_3d_geometric_check')
    components=[];models=[]
    for ci,a in enumerate(fields):
        checks=[];fcc=[];hpl=[];diamond=[]
        for level in LEVELS:
            labs,ids,top=topology(a>level);checks.append(dict(threshold=level,**top))
            f=fcc_model(labs,ids,A,a) if top['isolated'] else None;fcc.append(f)
            h=hpl_model(a,labs,ids,top,A,level);hpl.append(h)
            diamond.append(diamond_sphere_model(labs,ids,A) if top['isolated'] else None)
        components.append(dict(component='phi'+chr(65+ci),threshold_checks=checks))
        for name,tests in [('FCC',fcc),('HPL',hpl),('DSc',diamond)]:
            if all(t is not None for t in tests):models.append(dict(model=name,component='phi'+chr(65+ci),checks=tests))
    result=dict(candidate=None,components=components,models=models,reason='no_complete_model')
    if len(models)==1:
        m=models[0];candidate=m['model']+('i' if m['component']=='phiB' else '')
        # Third-component morphology needs joint interpretation; do not call C "inverse".
        if m['component'] in ('phiA','phiB'):result.update(candidate=candidate,reason='three_threshold_periodic_geometry_model',structural_component=m['component'])
    return result
def fuse(native,evidence):
    old=native.get('accepted_phase') or native.get('suggested_phase');candidate=evidence.get('candidate')
    if candidate:
        if candidate==old:return dict(status='source_and_geometry_supported',suggested_phase=old,accepted_phase=native.get('accepted_phase'),fusion='independent_geometry_agrees')
        return dict(status='review_geometry_supported_candidate',suggested_phase=candidate,accepted_phase=None,fusion='geometric_model_conflicts_with_or_extends_source')
    components=evidence.get('components',[])
    if old in {'DD','DG','SG'} and components and all(t['isolated'] for t in components[0]['threshold_checks']):
        return dict(status='review_structure_conflict',suggested_phase=None,accepted_phase=None,fusion='network_candidate_but_primary_component_isolated')
    spherical={'SC','SCi','BCC','BCCi','FCC','FCCi','A15','A15i','Sigma','Sigmai','DSc','DSci','HCP','iHPc','iHPa','iHPai','iHPs','PtS','PtSi'}
    if old in spherical and components and not any(all(t['isolated'] for t in c['threshold_checks']) for c in components):
        return dict(status='review_structure_conflict',suggested_phase=None,accepted_phase=None,fusion='sphere_candidate_without_persistent_isolated_domains')
    return dict(status=native['status'],suggested_phase=native.get('suggested_phase'),accepted_phase=native.get('accepted_phase'),fusion='source_retained_no_independent_model')
