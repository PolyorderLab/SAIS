"""Verified axis-supercell removal, before phase classification.
No phase names, cell-size prior or reference labels enter this transformation.
General sheared-supercell reduction is deliberately not implemented.
"""
import math
import numpy as np
VERSION='v4.0.2-experimental'
def reduce_supercell(fields,A,power_threshold=1e-12,relative_tolerance=1e-6):
    arrays=[np.asarray(x,dtype=float) for x in fields]
    shape=arrays[0].shape;d=len(shape)
    info=dict(version=VERSION,reduced=False,method='all-component verified axis translations',
              power_threshold=power_threshold,relative_tolerance=relative_tolerance,
              scope='integer repeats along saved cell vectors',axis_checks=[])
    if any(x.shape!=shape or not np.isfinite(x).all() for x in arrays):
        return fields,A,dict(info,reason='invalid_fields')
    transforms=[np.fft.fftn(x)/x.size for x in arrays]
    indices=[np.rint(np.fft.fftfreq(n)*n).astype(int) for n in shape]
    repeats=[]
    for ax in range(d):
        g=0
        for x,F in zip(arrays,transforms):
            p=abs(F.copy())**2;p[(0,)*d]=0
            spectrum=p.sum(axis=tuple(k for k in range(d) if k!=ax))
            peak=float(spectrum.max())
            if peak>1e-24:
                active=indices[ax][spectrum>power_threshold*peak]
                for h in active:g=math.gcd(g,abs(int(h)))
        g=max(1,g);errors=[]
        if g>1:
            reshape=[1]*d;reshape[ax]=shape[ax]
            phase=np.exp(2j*np.pi*indices[ax]/g).reshape(reshape)
            for x,F in zip(arrays,transforms):
                moved=np.fft.ifftn(F*phase).real*x.size
                errors.append(float(np.max(abs(moved-x)))/max(float(np.ptp(x)),1e-12))
        accepted=g>1 and max(errors,default=0)<=relative_tolerance
        info['axis_checks'].append(dict(axis=ax,proposed_repeats=g,translation_errors=errors,accepted=accepted))
        repeats.append(g if accepted else 1)
    info['repeats']=repeats
    if max(repeats)==1:return fields,A,dict(info,reason='no_verified_axis_supercell')
    grids=np.meshgrid(*indices,indexing='ij')
    keep=np.ones(shape,dtype=bool)
    for h,r in zip(grids,repeats):keep &= h%r==0
    errors=[]
    for x,F in zip(arrays,transforms):
        reconstructed=np.fft.ifftn(np.where(keep,F,0)).real*x.size
        errors.append(float(np.max(abs(reconstructed-x)))/max(float(np.ptp(x)),1e-12))
    info['reconstruction_errors']=errors
    if max(errors)>relative_tolerance:return fields,A,dict(info,reason='combined_translation_reconstruction_failed')
    # Keep every representable harmonic of the reduced period, including
    # negative Nyquist modes; odd output axes avoid merging +/- Nyquist bins.
    newshape=tuple(n if r==1 else 2*(n//2//r)+1 for n,r in zip(shape,repeats))
    mapped=tuple((h[keep]//r)%m for h,r,m in zip(grids,repeats,newshape))
    reduced=[]
    for F in transforms:
        small=np.zeros(newshape,dtype=complex)
        np.add.at(small,mapped,F[keep])
        reduced.append(np.fft.ifftn(small).real*np.prod(newshape))
    cell=np.asarray(A)@np.diag(1/np.array(repeats,dtype=float))
    info.update(reduced=True,reason='verified',original_shape=list(shape),reduced_shape=list(newshape),
                original_cell=np.asarray(A).tolist(),reduced_cell=cell.tolist())
    return reduced,cell,info

def select_axis_representation(fields,A,refs):
    """Keep the existing conventional cell unless a verified alternative improves
    the same scattering objective. No phase-specific preference is applied."""
    from scattering import identify
    reduced,cell,info=reduce_supercell(fields,A)
    if not info['reduced']:return fields,A,info
    if len(A) not in refs:return reduced,cell,info
    old=identify(fields[0],A,refs[len(A)])
    new=identify(reduced[0],cell,refs[len(A)])
    def score(result):
        return min((x['median_alignment'] for x in result.get('ranking',[])),default=float('inf'))
    oldscore,newscore=score(old),score(new)
    info.update(original_alignment=None if not np.isfinite(oldscore) else oldscore,
                reduced_alignment=None if not np.isfinite(newscore) else newscore)
    if newscore < oldscore-1e-10:
        info['selection']='verified_reduction_improves_reference_alignment'
        return reduced,cell,info
    info.update(reduced=False,verified_reduction_available=True,
                selection='retain_original_conventional_cell',reason='no_alignment_improvement')
    return fields,A,info

def select_representation(fields,A,refs):
    from scattering import identify,system
    from mixed_basis402 import candidate_cells,transform_fields,metric_compatible_references
    chosen,cell,info=select_axis_representation(fields,A,refs)
    if len(A) not in refs:return chosen,cell,info
    baseline=identify(chosen[0],cell,refs[len(A)])
    def score(r):
        return min((x['median_alignment'] for x in r.get('ranking',[])),default=float('inf'))
    bestscore=score(baseline)
    # Existing conventional metrics are already covered by axis reduction.
    # Mixed-basis search handles representations outside those metric families.
    from geometry061 import planar
    if bestscore<=1e-10 or (len(cell)==2 and planar(chosen,cell)):
        return chosen,cell,info
    if len(cell)==3:
        from structure_checks import topology
        layered=any(all(len(topology(v>t)[2]['spanning_axes'])==2 for t in (.35,.5,.65)) for v in chosen)
        if layered:
            from detail import describe
            evidence=describe(chosen,cell)
            if any(all(row.get('planes') and all(plane['isolated_holes'] and plane['hole_neighbors'] and all(n==4 for n in plane['hole_neighbors']) for plane in row['planes']) for row in component['checks']) for component in evidence):
                info['selection']='retain_verified_layer_geometry'
                return chosen,cell,info
    trials=[];best=None;best_names=set()
    for index,(multiplicity,kind,B,T) in enumerate(candidate_cells(fields,A)):
        values,check=transform_fields(fields,T)
        trial=dict(index=index,primitive_multiplicity=multiplicity,cell_system=kind,validation=check)
        if values is None:trials.append(trial);continue
        compatible=metric_compatible_references(B,refs[len(A)])
        trial['metric_compatible_templates']=list(compatible)
        if not compatible:
            trial['reason']='no_metric_compatible_reference';trials.append(trial);continue
        result=identify(values[0],B,compatible)
        value=score(result)
        trial.update(candidate=result.get('candidate'),alignment=None if not np.isfinite(value) else value)
        trials.append(trial)
        if value<bestscore-1e-10:
            bestscore=value;best=(values,B,T,check,index);best_names={result.get('candidate')}
        elif best is not None and abs(value-bestscore)<=1e-10:
            best_names.add(result.get('candidate'))
    info['mixed_basis_trials']=trials
    from morphology import canonical
    if len({canonical(n) for n in best_names})>1:
        info['mixed_basis_search']='conflicting_equal_score_candidates'
        return chosen,cell,info
    if best is not None:
        values,B,T,check,index=best
        info.update(reduced=True,selection='verified_integer_mixed_basis',
                    scope='axis repeats and bounded integer mixed bases',
                    selected_trial=index,basis_transform=T.tolist(),
                    reconstruction_errors=check['errors'],reduced_cell=B.tolist(),
                    reduced_shape=check['shape'])
        return values,B,info
    info['mixed_basis_search']='no_improving_validated_representation'
    return chosen,cell,info
