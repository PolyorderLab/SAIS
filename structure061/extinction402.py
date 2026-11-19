"""Diagnose finite-grid extinction conflicts without changing phase decisions.

For odd sample counts, h and h+N have opposite parity. Sampling can therefore
alias an allowed continuum reflection onto a forbidden discrete FFT index.
Neither interpolation nor projection onto allowed indices establishes the
unsampled continuum density. These diagnostics never waive forbidden power.
"""
import numpy as np

def diagnose(fields,phase,tolerance=1e-8):
 if phase not in ('BCC','FCC'):return None
 a=[np.asarray(v) for v in fields]
 if not a or any(v.ndim!=3 or not np.isfinite(v).all() for v in a):return None
 shape=np.array(a[0].shape)
 h=np.stack(np.meshgrid(*[np.rint(np.fft.fftfreq(n)*n).astype(int) for n in shape],indexing='ij'),axis=-1)
 allowed=(h.sum(axis=-1)%2==0) if phase=='BCC' else ((h[...,0]%2==h[...,1]%2)&(h[...,1]%2==h[...,2]%2))
 normq=np.max(abs(h)/(shape/2),axis=-1);q2=(h*h).sum(axis=-1)
 components=[]
 shifts=[[.5,.5,.5]] if phase=='BCC' else [[0,.5,.5],[.5,0,.5],[.5,.5,0]]
 for v in a:
  F=np.fft.fftn(v-v.mean())/v.size;p=abs(F)**2;total=float(p.sum())
  if total<=1e-30:return None
  bad=p*(~allowed);bp=float(bad.sum());sh=[]
  for t in shifts:
   error=np.fft.ifftn(F*(np.exp(2j*np.pi*(h*np.array(t)).sum(axis=-1))-1))*v.size
   sh.append(dict(translation=t,rms=float(np.sqrt(np.mean(abs(error)**2))),
    method='exact_grid_shift' if np.allclose(shape*np.array(t),np.rint(shape*np.array(t))) else 'bandlimited_interpolant_not_independent_measurement'))
  indices=np.argsort(bad.ravel())[-4:][::-1]
  components.append(dict(forbidden_fraction=bp/total,
   low_q_fraction=float(bad[q2<=16].sum()/total),
   interior_fraction=float(bad[normq<=.5].sum()/total),
   outer_band_share=float(bad[normq>.5].sum()/bp) if bp else 0.,
   strongest_modes=[dict(h=h.reshape(-1,3)[i].tolist(),fraction=float(bad.ravel()[i]/total)) for i in indices],
   translation_mismatch=sh))
 fraction=max(c['forbidden_fraction'] for c in components)
 odd=bool(np.any(shape%2))
 tail=all(c['outer_band_share']>=.9 for c in components)
 if fraction<=tolerance:category='centering_within_strict_gate'
 elif odd:category='odd_grid_tail_dominated_ambiguity' if tail else 'odd_grid_low_and_high_frequency_ambiguity'
 else:category='commensurate_grid_symmetry_violation_origin_unresolved'
 return dict(phase_tested=phase,grid_shape=shape.tolist(),half_cell_grid_commensurate=not odd,
  category=category,components=components,
  physical_distortion_proven=False,numerical_origin_proven=False,
  decision_effect='diagnostic_only_no_threshold_waiver',
  next_validation='same_cell_same_system_independent_even_finer_grid_solve_and_tighter_residual',
  caveat='A single sampled field cannot establish continuum symmetry or thermodynamic stability. Low-q leakage can also be aliased; residual tolerance is not a density-error bound.')

def annotate(result,fields):
 ev=result.get('particle_evidence')
 if not ev:return result
 phases=[n for n in ev.get('arrangements',[]) if n in ('BCC','FCC')]
 if not phases:return result
 out=dict(result);out['extinction_diagnostics']=[diagnose(fields,n,ev['policy']['centering_power_tolerance']) for n in phases]
 return out
