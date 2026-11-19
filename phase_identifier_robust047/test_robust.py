"""Real-field controls for the new basis and weak-peak evidence."""
from pathlib import Path
import unittest
import numpy as np

from classify_batch import load
from robust_checks import hex_orthogonal_basis, spectral_weak_peak_audit

ROOT = Path(__file__).resolve().parents[1]


class RealFieldControls(unittest.TestCase):
    def test_hex_basis_requires_primitive_translation(self):
        fields, A, _ = load(ROOT/'fcc_review033/inputs/run2720')
        result = hex_orthogonal_basis(A, 2, fields)
        self.assertEqual(result['determinant'], 2)
        self.assertLess(max(result['component_half_translation_normalized_rms']), .002)
        bad = [f.copy() for f in fields]
        x = np.arange(bad[0].shape[0])[:, None, None]
        bad[0] += .1*np.cos(2*np.pi*x/bad[0].shape[0])
        self.assertIsNone(hex_orthogonal_basis(A, 2, bad))

    def test_weak_fcc_peak_differs_from_dd(self):
        f20, a20, _ = load(ROOT/'fcc_review033/inputs/run8442')
        result20 = spectral_weak_peak_audit(f20[0], a20)
        self.assertGreater(result20['relative_threshold_sweep'][0]['forbidden_above'], 0)
        self.assertEqual(result20['relative_threshold_sweep'][1]['forbidden_above'], 0)
        fdd, add, _ = load(ROOT/'catalog_benchmark037/inputs/DD.mat')
        resultdd = spectral_weak_peak_audit(fdd[0], add)
        self.assertGreater(resultdd['relative_threshold_sweep'][1]['forbidden_above'], 0)


if __name__ == '__main__':
    unittest.main()
