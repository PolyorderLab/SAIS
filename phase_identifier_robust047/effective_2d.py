"""Conservative all-component dimensional routing, independent of phase names."""
from dataclasses import dataclass,asdict
import numpy as np

@dataclass(frozen=True)
class Config:
    absolute_tolerance:float=1e-10
    relative_tolerance:float=1e-10
    geometry_tolerance:float=1e-10

def project(fields,A,config=Config()):
    A=np.asarray(A,dtype=float)
    if not np.isfinite(A).all() or A.ndim!=2 or A.shape[0]!=A.shape[1] or abs(np.linalg.det(A))<1e-14:
        raise ValueError('Invalid cell basis')
    if any(a.ndim!=A.shape[0] or a.shape!=fields[0].shape or not np.isfinite(a).all() for a in fields):raise ValueError('Invalid fields')
    info={'original_dimension':len(A),'effective_dimension':len(A),'reduced':False,'config':asdict(config)}
    if len(A)!=3:return fields,A,dict(info,reason='input_not_3d')
    spans=[float(np.ptp(a)) for a in fields]
    variations=[[float(np.max(np.ptp(a,axis=axis))) for a in fields] for axis in range(3)]
    limits=[config.absolute_tolerance+config.relative_tolerance*s for s in spans]
    axes=[axis for axis in range(3) if fields[0].shape[axis]>=2 and all(v<=t for v,t in zip(variations[axis],limits))]
    info.update(axis_variations=variations,component_limits=limits,invariant_axes=axes)
    if len(axes)!=1:return fields,A,dict(info,reason='requires_exactly_one_invariant_axis')
    axis=axes[0];keep=[i for i in range(3) if i!=axis]
    lengths=np.linalg.norm(A,axis=0)
    cross=max(abs(float(A[:,axis]@A[:,j]))/(lengths[axis]*lengths[j]) for j in keep)
    info['axis_plane_cosine_max']=cross
    if cross>config.geometry_tolerance:return fields,A,dict(info,reason='oblique_axis_not_supported')
    basis=A[:,keep];gram=basis.T@basis;A2=np.linalg.cholesky(gram).T
    reduced=[a.mean(axis=axis) for a in fields]
    assert np.allclose(A2.T@A2,gram,rtol=1e-12,atol=1e-12)
    errors=[float(np.max(np.abs(a-np.expand_dims(b,axis=axis)))) for a,b in zip(fields,reduced)]
    info.update(reduced=True,effective_dimension=2,removed_axis=axis,retained_axes=keep,
                original_plane_gram=gram.tolist(),projected_basis=A2.tolist(),reconstruction_errors=errors,
                method='all-component axial mean; no smoothing in retained plane',reason='resolved_2d_in_orthogonal_extrusion')
    return reduced,A2,info
