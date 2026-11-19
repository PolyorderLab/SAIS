"""Experimental all-component amplitude guard; never modifies input arrays."""
from dataclasses import dataclass,asdict
import numpy as np

@dataclass(frozen=True)
class GuardConfig:
    dis_ptp: float=1e-5
    dis_std: float=1e-6
    review_ptp: float=1e-4
    review_std: float=1e-5
    fraction_sum_tolerance: float=1e-3

def evaluate(fields,config=GuardConfig()):
    if not 0<config.dis_ptp<config.review_ptp or not 0<config.dis_std<config.review_std:
        raise ValueError('Invalid guard thresholds')
    arrays=[np.asarray(x,dtype=np.float64) for x in fields]
    result={'version':'all-component-amplitude-guard-v1','config':asdict(config),
            'thresholds_are_validated_physical_phase_boundary':False}
    if len(arrays) not in (2,3) or not arrays[0].size or len({a.shape for a in arrays})!=1 or not all(np.isfinite(a).all() for a in arrays):
        return dict(result,decision='review_invalid_input',reason='invalid component count, shape or nonfinite values')
    stats=[dict(component='phi'+chr(65+i),min=float(a.min()),max=float(a.max()),mean=float(a.mean()),
                ptp=float(np.ptp(a)),std=float(a.std())) for i,a in enumerate(arrays)]
    error=float(np.max(np.abs(sum(arrays)-1)))
    result.update(components=stats,total_component_fraction_error=error,dimension=arrays[0].ndim,
                  peak_to_peak_max=max(s['ptp'] for s in stats),std_max=max(s['std'] for s in stats))
    if error>config.fraction_sum_tolerance:
        return dict(result,decision='review_input_quality',reason='component fractions do not sum to one within tolerance')
    if all(s['ptp']<=config.dis_ptp and s['std']<=config.dis_std for s in stats):
        return dict(result,decision='DIS_candidate',reason='all components are within the declared numerical amplitude limits')
    if all(s['ptp']<=config.review_ptp and s['std']<=config.review_std for s in stats):
        return dict(result,decision='review_weak_amplitude',reason='weak-amplitude boundary region; do not force DIS or ordered label')
    return dict(result,decision='native_required',reason='resolved density variation; retain original identification workflow')
