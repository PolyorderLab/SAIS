"""Conservative evidence gate. Proposals only; never writes training/database labels.

Energy convention: JuPolyorder F_DIS = enthalpy_FH + F_ig.
Only single binary BCP melts are supported here. Review bands are experimental,
not a physical ODT criterion or a rigorous estimate of numerical energy error.
"""
from dataclasses import dataclass, asdict
from pathlib import Path
import hashlib
import csv
import math
import numpy as np
import yaml

@dataclass(frozen=True)
class Policy:
    dis_ptp: float = 1e-5
    dis_std: float = 1e-6
    weak_ptp: float = 1e-2
    weak_std: float = 3e-3
    energy_review_tolerance: float = 1e-6
    composition_tolerance: float = 1e-4
    density_tolerance: float = 1e-3

def read_yaml(path, hashes):
    if not path.is_file(): return None
    data = path.read_bytes()
    hashes[str(path)] = hashlib.sha256(data).hexdigest()
    value = yaml.safe_load(data.decode('utf-8-sig'))
    if not isinstance(value, dict): raise ValueError('YAML must contain a mapping')
    return value

def finite_number(value):
    x = float(value)
    if not math.isfinite(x): raise ValueError('nonfinite scalar')
    return x

def dis_baseline(config, fields, policy):
    system = config['system']
    components = system['components']
    if len(components) != 1 or len(fields) != 2:
        raise ValueError('unsupported mixture or component count')
    c = components[0]
    if c['type'] != 'BCP': raise ValueError('unsupported molecule type')
    vf, alpha, density = map(finite_number, (c['volume_fraction'], c['length'], system['chain_density']))
    if abs(vf-1) > 1e-10 or alpha <= 0 or density <= 0:
        raise ValueError('not a supported single-polymer melt')
    fractions = {}
    for b in c['blocks']:
        length = finite_number(b['length'])
        if length <= 0: raise ValueError('nonpositive block length')
        fractions[b['segment']] = fractions.get(b['segment'], 0.) + length
    if set(fractions) != {'A', 'B'} or abs(sum(fractions.values())-1) > 1e-8:
        raise ValueError('unsupported species or block-length normalization')
    pairs = system['χN_map']
    if len(pairs) != 1 or set(pairs[0][:2]) != {'A', 'B'}:
        raise ValueError('unsupported interaction map')
    chi = finite_number(pairs[0][2])
    if any(abs(float(a.mean())-fractions[s]) > policy.composition_tolerance
           for a, s in zip(fields, ('A', 'B'))):
        raise ValueError('configuration composition does not match density')
    ideal = vf/alpha * (math.log(density*vf/alpha)-1)
    return dict(value=chi*fractions['A']*fractions['B']+ideal,
                convention='JuPolyorder F_DIS = enthalpy_FH + F_ig',
                f_A=fractions['A'], f_B=fractions['B'], chi_N=chi,
                ideal_gas_term=ideal, relative_chain_length=alpha,
                chain_density=density)

def convergence_evidence(summary, config):
    if summary is None or config is None:
        return dict(passed=False, reason='missing_summary_or_config')
    try:
        s, c = config['scft'], config['cellopt']
        residual = finite_number(summary['residual'])
        stress = finite_number(summary['stress_norm'])
        tol, gtol = finite_number(s['tol']), finite_number(c['gtol'])
        # Only the residual stopping convention has been checked against sources.
        mode_ok = s.get('tolmode') == 'Residual' and s.get('relative') is False
        status_ok = summary.get('convergence') == 'JuPolyorder.Successful()'
        valid = residual >= 0 and stress >= 0 and tol > 0 and gtol > 0
        passed = valid and mode_ok and status_ok and residual <= tol and stress <= gtol
        return dict(passed=bool(passed), reported_status=summary.get('convergence'),
                    residual=residual, residual_tolerance=tol, stress_norm=stress,
                    stress_tolerance=gtol, supported_stopping_convention=mode_ok,
                    reason='summary_and_thresholds_pass' if passed else 'convergence_not_verified')
    except (ValueError, TypeError, KeyError):
        return dict(passed=False, reason='invalid_or_missing_convergence_metadata')

def trace_evidence(metadata, hashes):
    result = {'amplitude_history_available': False}
    for name in ('trace_cell.csv', 'trace_solve.csv'):
        path = metadata/name
        if not path.is_file(): continue
        data = path.read_bytes(); hashes[str(path)] = hashlib.sha256(data).hexdigest()
        try:
            lines = data.decode('utf-8-sig').splitlines()
            rows = list(csv.DictReader(lines, delimiter='\t'))
            last = rows[-1] if rows else {}
            result[name] = dict(records=len(rows), last_record=last,
                                usage='context_only_not_final_amplitude_or_stopping_proof')
        except (UnicodeError, csv.Error):
            result[name] = {'reason': 'unreadable_trace'}
    return result

def apply(result, fields, source, metadata_dir=None, policy=Policy()):
    """Return a new result without mutating result, density arrays, or source files."""
    if not (0 < policy.dis_ptp < policy.weak_ptp and
            0 < policy.dis_std < policy.weak_std and policy.energy_review_tolerance > 0):
        raise ValueError('Invalid review policy')
    source = Path(source); metadata = Path(metadata_dir) if metadata_dir else source
    hashes = {}; errors = []; summary = config = None
    try:
        if source.is_dir(): summary = read_yaml(source/'summary.yml', hashes)
        if metadata.is_dir(): config = read_yaml(metadata/'config.yml', hashes)
    except (ValueError, UnicodeError, yaml.YAMLError) as exc:
        errors.append(str(exc))
    a = [np.asarray(x, dtype=float) for x in fields]
    valid = (len(a) == 2 and all(x.size and np.isfinite(x).all() for x in a)
             and len({x.shape for x in a}) == 1)
    stats = [dict(min=float(x.min()), max=float(x.max()), mean=float(x.mean()),
                  ptp=float(np.ptp(x)), std=float(x.std())) for x in a] if valid else []
    if valid:
        valid = (np.max(abs(sum(a)-1)) <= policy.density_tolerance and
                 all(x.min() >= -policy.density_tolerance and x.max() <= 1+policy.density_tolerance for x in a))
    uniform = bool(valid and all(x['ptp'] <= policy.dis_ptp and x['std'] <= policy.dis_std for x in stats))
    weak = bool(valid and all(x['ptp'] <= policy.weak_ptp and x['std'] <= policy.weak_std for x in stats))
    conv = convergence_evidence(summary, config)
    energy = {'available': False}
    if valid and summary is not None and config is not None:
        try:
            base = dis_baseline(config, a, policy)
            value = finite_number(summary['F'])
            energy = dict(available=True, F=value, baseline=base, delta_F=value-base['value'],
                          tolerance=policy.energy_review_tolerance,
                          tolerance_is_physical_or_rigorous_error_bound=False)
        except (KeyError, ValueError, TypeError) as exc:
            energy['reason'] = str(exc)
    else: energy['reason'] = 'missing_metadata_or_invalid_density'
    morphology = result.get('candidate')
    if result.get('status') == 'DIS_candidate': morphology = 'DIS'
    proposal = morphology
    if not valid: status = 'review_invalid_density'
    elif config is None or summary is None or not energy['available']:
        status = 'review_missing_or_unsupported_metadata'
    elif not conv['passed']: status = 'review_convergence_unverified'
    elif uniform:
        status = 'joint_DIS_supported' if abs(energy['delta_F']) <= policy.energy_review_tolerance else 'review_uniform_energy_conflict'
        if status == 'joint_DIS_supported': proposal = 'DIS'
    elif weak and energy['delta_F'] >= -policy.energy_review_tolerance:
        status = 'review_DIS_vs_weak_order'
    elif energy['delta_F'] > policy.energy_review_tolerance:
        status = 'review_above_DIS_energy'
    elif not morphology or str(result.get('status', '')).startswith('review_'):
        status = 'review_structure_unresolved'
    elif weak: status = 'review_weak_order_below_DIS'
    elif energy['delta_F'] >= -policy.energy_review_tolerance:
        status = 'review_ordered_energy_near_DIS'
    else: status = 'joint_ordered_supported'
    supported = status in ('joint_DIS_supported', 'joint_ordered_supported')
    out = dict(result)
    out.update(version='061-joint-evidence-experimental', morphology_candidate=morphology,
               morphology_status=result.get('status'), candidate=proposal if supported else None,
               physical_phase_candidate=morphology, status=status,
               review_category=('numerical_quality' if not valid else 'convergence' if not conv['passed'] else 'near_DIS' if weak else 'stability' if status=='review_above_DIS_energy' else 'structure' if not supported else None),
               confirmed_label=None, database_written=False, training_label_written=False)
    out['joint_evidence'] = dict(policy=asdict(policy), density_valid=bool(valid),
        amplitude=stats, numerically_uniform=uniform, weak_amplitude=weak,
        convergence=conv, energy=energy, metadata_errors=errors,
        traces=trace_evidence(metadata, hashes) if metadata.is_dir() else {},
        metadata_hashes=hashes, thermodynamic_global_stability_proven=False,
        interpretation='Proposal only. Lower than DIS does not establish global stability.')
    return out
