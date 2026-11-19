"""Basis and normalized reciprocal-space checks for SCFT fields."""
import numpy as np

def hex_orthogonal_basis(A: np.ndarray, normal: int, fields):
    """Find the conventional index-two orthogonal cell of a 120-degree hex net."""
    if A.shape != (3, 3):
        return None
    axes = [i for i in range(3) if i != normal]
    if any(abs(np.dot(A[:, normal], A[:, j])) /
           (np.linalg.norm(A[:, normal]) * np.linalg.norm(A[:, j])) > 1e-4
           for j in axes):
        return None
    for short, long in (axes, axes[::-1]):
        u, v = A[:, short], A[:, long]
        ratio = np.linalg.norm(v) / np.linalg.norm(u)
        cosine = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
        if abs(ratio - np.sqrt(3)) > 0.01 or abs(cosine) > 1e-4:
            continue
        h1, h2 = u, (v-u)/2
        cos_hex = np.dot(h1, h2)/(np.linalg.norm(h1)*np.linalg.norm(h2))
        if abs(cos_hex+0.5) > 0.01:
            continue
        if any(f.shape[short] % 2 or f.shape[long] % 2 for f in fields):
            continue
        errors = [float(np.sqrt(np.mean((f-np.roll(f,
                    (f.shape[short]//2, f.shape[long]//2),
                    (short, long)))**2))/max(float(np.std(f)), 1e-12))
                  for f in fields]
        if max(errors) > 0.02:
            continue
        return dict(status="supported_index_two_hex_basis", short_axis=short,
                    long_axis=long, layer_normal_axis=normal,
                    orthogonal_to_hex_columns=[[1, 1], [0, 2]],
                    determinant=2, orthogonal_edge_ratio=float(ratio),
                    derived_hex_cosine=float(cos_hex),
                    component_half_translation_normalized_rms=errors,
                    note="metric and primitive-translation evidence; morphology is checked separately")
    return None


def spectral_weak_peak_audit(field: np.ndarray, A: np.ndarray):
    """Normalize FFT power and expose low-order forbidden-peak sensitivity."""
    if field.ndim != 3 or A.shape != (3, 3):
        return dict(status="not_applicable")
    g = A.T @ A
    edge = np.sqrt(np.diag(g))
    if not np.allclose(g/np.outer(edge, edge), np.eye(3), atol=1e-4) or \
            np.ptp(edge)/edge.mean() > 1e-3:
        return dict(status="not_cubic_conventional_cell")
    power = np.abs(np.fft.fftn(field-field.mean())/field.size)**2
    h = [np.rint(np.fft.fftfreq(n)*n).astype(int) for n in field.shape]
    H = np.meshgrid(*h, indexing="ij")
    low = np.logical_and.reduce([np.abs(x) <= 4 for x in H])
    nonzero = np.logical_or.reduce([x != 0 for x in H])
    parity = [x % 2 for x in H]
    fcc_allowed = (parity[0] == parity[1]) & (parity[1] == parity[2])
    usable = low & nonzero
    peak = float(power[usable].max(initial=0))
    forbidden = power[usable & ~fcc_allowed]
    allowed = power[usable & fcc_allowed]
    if peak <= 0:
        return dict(status="no_low_order_power")
    ratios = [1e-6, 1e-5, 1e-4]
    return dict(status="evaluated", normalized_peak_power=peak,
                max_forbidden_to_main=float(forbidden.max(initial=0)/peak),
                total_forbidden_to_allowed=float(forbidden.sum()/max(allowed.sum(), 1e-300)),
                relative_threshold_sweep=[dict(threshold=t,
                    forbidden_above=int(np.count_nonzero(forbidden > t*peak)),
                    allowed_above=int(np.count_nonzero(allowed > t*peak))) for t in ratios],
                threshold_near_peak=bool(np.any((forbidden/peak > 1e-6) &
                                                (forbidden/peak < 1e-4))),
                note="normalized FFT power; this check is FCC-specific")


