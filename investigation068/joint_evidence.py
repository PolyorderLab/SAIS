"""068: require evaluation-index alignment before comparing trace residuals.

Trace samples may precede final refinement. A conflict requests verification,
not a claim that the original solver failed. Missing traces are explicit.
"""
import ast
import math
from base061 import apply as apply061

def vector(value):
    try:
        x = ast.literal_eval(value) if isinstance(value, str) else value
        if not isinstance(x, (list, tuple)) or not x:
            return None
        x = [float(v) for v in x]
        return x if all(math.isfinite(v) for v in x) else None
    except (ValueError, TypeError, SyntaxError):
        return None

def audit(traces, convergence):
    cell = traces.get('trace_cell.csv', {}).get('last_record', {})
    solve = traces.get('trace_solve.csv', {}).get('last_record', {})
    a, b = vector(cell.get('unitcell')), vector(solve.get('unit cell'))
    result = dict(status='not_comparable', conflict=False,
                  interpretation='Trace samples are not guaranteed to be final solver state.')
    try:
        trace_eval = int(solve['nevals'])
        final_eval = int(cell['nevals_solve'])
    except (KeyError, ValueError, TypeError):
        result['status'] = 'trace_evaluation_alignment_unknown'
        return result
    result.update(trace_evaluation=trace_eval, final_evaluation=final_eval)
    if trace_eval != final_eval:
        result['status'] = 'trace_evaluation_mismatch_not_final_state'
        return result
    if a is None or b is None or len(a) != len(b):
        return result
    same = all(math.isclose(x, y, rel_tol=1e-9, abs_tol=1e-12) for x, y in zip(a, b))
    result['same_final_cell'] = same
    if not same:
        return result
    try:
        residual = float(solve['residual'])
        tol = float(convergence['residual_tolerance'])
        if not math.isfinite(residual) or residual < 0 or not math.isfinite(tol) or tol <= 0:
            raise ValueError()
    except (KeyError, ValueError, TypeError):
        result.update(status='invalid_trace_residual', conflict=True)
        return result
    result.update(trace_residual=residual, residual_tolerance=tol,
                  trace_residual_over_tolerance=residual/tol)
    conflict = bool(convergence.get('passed') and residual > tol)
    result.update(status='final_cell_trace_requires_verification' if conflict else 'no_residual_conflict_detected',
                  conflict=conflict)
    return result

def apply(*args, **kwargs):
    out = apply061(*args, **kwargs)
    e = out['joint_evidence']
    check = audit(e['traces'], e['convergence'])
    e['trace_consistency'] = check
    out['version'] = '068-aligned-trace-experimental'
    out['ruleset_revision'] = 'v4.0.3-polyorder-cell-scope-20261124'
    if check['conflict']:
        out['prior_gate_status'] = out['status']
        out['prior_gate_candidate'] = out['candidate']
        out.update(candidate=None, status='review_convergence_trace_conflict', review_category='convergence')
    return out
